#!/usr/bin/env python3
"""liftoff_cross — the cross-species Liftoff runs of Phase E, tier B (docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md section 7, steps 2-3).

The gorilla annotation is the REFERENCE and each outgroup assembly (H = CHM13v2.0, C = mPanTro3, O = mPonPyg2) the TARGET of a `liftoff -copies` run that holds
at most one record of every called pool component (the rank-sharded plan of expansion_table.rank_plan). Everything here is code and tests; NOTHING in this
file is executed by the repository's tests against Liftoff or minimap2, and the `run` command refuses to start without the commit hash of the Amendment 2
that authorises it (a heavy 'go' of the user is a precondition of that Amendment).

Gate 7 comes first (section 7): the two controls are lifted and checked BEFORE any verdict-set record is lifted, in calls of their own. Three groups of runs,
in this order, and a call never mixes groups: `pos` = the rank plan of the LRPAP1 component and the fragment component (labels pos1..), `neg` = the 100
planted records in two runs (ctl1, ctl2), `verdict` = the rank plan of the other called components (r01..). The gate is asked before the first unit of a group
starts: `neg` needs the positive control to pass, `verdict` needs both. The verdict set holds no pool component of the controls.

Reused from figures/_liftoff.py (imported, not copied): `parse_lifted` (a Liftoff GFF3 -> parent-level loci.tsv rows), `mm2_wrapper` (the per-call minimap2
result cache Liftoff's `-m` takes), `TIME_FMT`, and the registered parameters LIFT_A / LIFT_S / LIFT_SC / LIFT_OVERLAP / LIFT_THREADS / LIFT_TYPES /
HARD_LIMIT_S / LIFTOFF_ENV. New here: the cross-species target (the reference FASTA differs from the target FASTA), the one-pass GFF router for arbitrary run
labels, `timed_watch` (`_timed` plus a memory watchdog that stops a run while it runs), the budgeted driver with its call deadline and stop rules, the Gate 7
sequencing, and the build of the E-table from the finished runs.

Stop rules (all end the call, none is retried; the state on disk stays valid): a run killed at the hard limit is split at the record-free cut nearest its middle and
its halves run in the next call (no cut: stop); a summed RSS above 16 GB (watchdog); a non-zero Liftoff exit; a run whose parent records are not all lifted or listed
unmapped (Liftoff does not look at the exit status of its minimap2); the summed Liftoff wall of one outgroup above twice its predicted wall (stop rule 6, E_ROOT/wall.json);
Gate 7 failing before the `neg` or `verdict` group (the controls are evaluated every time); a liftoff / minimap2 of this root still alive, or another driver (or the leftovers
of a killed one) holding E_ROOT/run.lock; a plan, code, tool, command or input other than the one the first call froze (E_ROOT/plan.digest, run_identity.json); a control or
verdict run that is not finished, or was made by another command or from another in.gff, when the table is built.

The Liftoff tree never outlives the driver: it runs under `timeout -s KILL` inside its own session, a SIGTERM / SIGHUP / SIGINT to the driver kills the session first
(from /proc, without `ps`), the session is verified empty before every return, and the Liftoff main process and its pool workers inherit the driver's run.lock (minimap2
and the cache wrapper do not; the stray guard covers them). Nothing is frozen before the first Liftoff process is about to start.

CLI:  python3 -B bench/hierarchy/liftoff_cross.py plan   [--hub-mode nodes|records]       (light: the plans, the schedule and the wall budget, no Liftoff)
      python3 -B bench/hierarchy/liftoff_cross.py gffs   [--hub-mode ...] [--out DIR]     (light, about 10 s: writes the per-run reference GFFs and their assignment marker)
      python3 -B bench/hierarchy/liftoff_cross.py run --amendment-commit HASH --budget-s 530 [--hub-mode ...]   (HEAVY: RLOCK_WAIT=5 RLOCK_TIMEOUT=585 tools/rlock.sh heavy, one call)
      python3 -B bench/hierarchy/liftoff_cross.py build [--out DIR] [--table-dir DIR] [--strict-single-copy] [--rehearsal]  (light: Gate 7 from the finished control runs, then the E-table)
Exit codes: 0 done, 75 call again, 3 a stop rule (STOP: ...), 2 refused (no commit hash, a rehearsal flag with the real --out, a build into a non-empty directory, ...),
4 build wrote an INVALID table, 70 a defect of this code (traceback), 128 + n a signal. The options that change what a run is made of (--liftoff-bin, --minimap2,
--no-code-check, --lift-root, --gorilla-gff, --gorilla-fasta, --n-planted, run --strict-single-copy, build --rehearsal) are for rehearsals: refused when --out is the real E_ROOT.
"""
from __future__ import annotations

import argparse
import collections
import contextlib
import csv
import fcntl
import gzip
import hashlib
import json
import os
import re
import shutil
import signal
import struct
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path

CALL_T0 = time.monotonic()          # the start of this process: every call limit below is counted from here, so the start-up work (plan, hashing, GFFs) is charged
HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "figures"))
sys.path.insert(0, str(HERE))
import _liftoff as L  # noqa: E402  (figures/_liftoff.py: the self-lift code, reused)
import expansion_table as X  # noqa: E402

E_ROOT = "/mnt/linuxdisk/tmp/hier_run/e"
HOME = "/mnt/linuxdisk/home/juanfraitu"
GORILLA_FASTA = HOME + "/_from_wsl/winloci_scratch/GGO.fasta"
GORILLA_GFF = HOME + "/winloci_data/GGO_genomic.gff"
LIFTOFF_BIN = L.LIFTOFF_ENV + "/bin/liftoff"
MINIMAP2_ENV = L.LIFTOFF_ENV + "/bin/minimap2"

LIFT_A, LIFT_S, LIFT_SC, LIFT_OVERLAP, LIFT_THREADS = L.LIFT_A, L.LIFT_S, L.LIFT_SC, L.LIFT_OVERLAP, L.LIFT_THREADS
HARD_LIMIT_S = L.HARD_LIMIT_S                      # 590: a Liftoff call is killed after this (the 10-minute rule)
CALL_BUDGET_S = 530                                # a call starts no unit whose predicted time (counted from the start of the PROCESS) does not fit in this
TERM_GRACE_S = 10.0                                # SIGTERM to the Liftoff session, then SIGKILL after this many seconds
KILL_WAIT_S = 10.0                                 # after SIGKILL, wait at most this long for the session to be empty (a process in uninterruptible I/O can take a while)
CALL_LIMIT_S = 555                                 # no Liftoff process outlives this many seconds from the START OF THE PROCESS: 555 + the worst kill sequence (TERM_GRACE_S +
#                                                    KILL_WAIT_S = 20 s) + the clean-up of inter/ and the bookkeeping (about 5 s) < 585 s = RLOCK_TIMEOUT, which is < 590 s
INSESSION_MARGIN_S = 5.0                           # the `timeout -s KILL` inside the Liftoff session fires this long after the driver's own hard limit (a backstop for a dead driver)
POLL_S = 0.5                                       # the watchdog polls the summed RSS of the session this often
WALL_FACTOR = 2.0                                  # stop rule 6: the summed Liftoff wall of one outgroup may not pass this multiple of its predicted wall
TIME_BIN = "/usr/bin/time"
TIMEOUT_BIN = "/usr/bin/timeout"
N_PLANTED = 100                                    # the planted negative control (registered: 100 records in 50 fake pairs); only the tests change it
INPUT_BIG_BYTES = 10 ** 9                          # an input below this size is pinned by its sha256 as well, above it by size and mtime (its sha256 is Gate 0's)
SHALLOW_SHA_BYTES = 100 * 10 ** 6                  # the calls after the first re-hash only inputs below this (the 693 MB annotation is hashed by the first call and by build)
AMENDMENT_FILE = "docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md"
AMENDMENT_MARK = "## Amendment 2"
REHEARSAL_MARKER = "REHEARSAL.json"               # written into --out by a run that used a stub flag: `build` refuses such a directory unless told it is a rehearsal
RSS_LIMIT_GB = 16.0                                # stop rule: a run above this (the self-lifts peaked at 10.2 to 13.8 GB)
PER_RECORD_S = 0.09                                # upper end of the measured 0.05 to 0.09 s per record (c12_shardfit.out)
STARTUP_S = 25.0                                   # index load and set-up of one run (10 to 20 s measured, plus margin)

OUTGROUPS = {
    "H": dict(species="human", assembly="CHM13v2.0 (RefSeq RS_2025_08)", genome=HOME + "/winloci_data/chm13v2.0.fa",
              mmi=HOME + "/winloci_data/chm13v2.0.fa.mmi", selflift=X.LIFT_ROOT + "/human/liftoff_loci.tsv"),
    "C": dict(species="chimpanzee", assembly="mPanTro3 (primary)",
              genome=HOME + "/winloci_data/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna", mmi=HOME + "/ggo_npip/ref/PTR.fa.mmi",
              selflift=X.LIFT_ROOT + "/chimpanzee/liftoff_loci.tsv"),
    "O": dict(species="orangutan", assembly="mPonPyg2 (primary)",
              genome=HOME + "/winloci_data/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna",
              mmi=HOME + "/o1_falsemerge/lit/evid/PPY/genome/PPY.fa.mmi", selflift=X.LIFT_ROOT + "/orangutan/liftoff_loci.tsv"),
}


class StopRun(RuntimeError):
    """a stop rule fired: do not retry, report to the user"""


# ----------------------------------------------------------------------------------------------- the command
def types_text() -> str:
    """Liftoff lifts GFF features of type gene only unless -f names others: pseudogene records are lifted too"""
    return "\n".join(L.LIFT_TYPES) + "\n"


def liftoff_command(target: str, reference: str, minimap2: str, liftoff_bin: str = LIFTOFF_BIN, threads: int = LIFT_THREADS) -> list[str]:
    """the registered command (section 7, step 2), run in the run directory that holds in.gff and types.txt; TARGET is the outgroup, REFERENCE gorilla"""
    return [liftoff_bin, str(target), str(reference), "-g", "in.gff", "-f", "types.txt", "-o", "lifted.gff3", "-u", "unmapped.txt", "-dir", "inter",
            "-copies", "-sc", str(LIFT_SC), "-a", str(LIFT_A), "-s", str(LIFT_S), "-overlap", str(LIFT_OVERLAP), "-p", str(threads), "-m", str(minimap2)]


# ----------------------------------------------------------------------------------------------- the reference GFF of every run, in one pass
_ID = re.compile(r"(?:^|;)ID=([^;]+)")
_PARENT = re.compile(r"(?:^|;)Parent=([^;]+)")


def write_run_gffs(gff_path: str, assign: dict[str, str], out_dir: str) -> dict[str, int]:
    """Write `out_dir/<label>/in.gff`: the lines of every record assigned to <label> and of their descendants, in file order (RefSeq GFF3 lists parents
    before children; `region`, `cDNA_match` and `match` lines are left out, as figures/_liftoff.py::prepare does). One pass; the routing follows the
    FIRST parent like prepare(). Returns the number of records written per label; a record that is not in the annotation is a KeyError."""
    top_of: dict[str, str] = {}
    handles, written = {}, collections.defaultdict(set)
    try:
        for label in sorted(set(assign.values())):
            d = os.path.join(out_dir, label)
            os.makedirs(d, exist_ok=True)
            handles[label] = open(os.path.join(d, "in.gff.tmp"), "w")
            handles[label].write("##gff-version 3\n")
        with open(gff_path) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.split("\t", 8)
                if len(f) < 9 or f[2] in ("cDNA_match", "match"):
                    continue
                col9 = f[8].rstrip("\n")
                mp = _PARENT.search(col9)
                if mp is None:
                    mi = _ID.search(col9)
                    if mi is None:
                        continue
                    rid = mi.group(1)
                    if rid in assign:
                        top_of[rid] = rid
                        handles[assign[rid]].write(line)
                        written[assign[rid]].add(rid)
                    continue
                top = top_of.get(mp.group(1).split(",")[0])
                if top is None:
                    continue
                mi = _ID.search(col9)
                if mi is not None:
                    top_of[mi.group(1)] = top
                handles[assign[top]].write(line)
    finally:
        for h in handles.values():
            h.close()
    missing = sorted(set(assign) - set(top_of))
    if missing:
        raise KeyError(f"{len(missing)} record(s) not in {gff_path}: {missing[:5]}")
    for label in handles:
        d = os.path.join(out_dir, label)
        os.replace(os.path.join(d, "in.gff.tmp"), os.path.join(d, "in.gff"))
    return {label: len(written[label]) for label in sorted(handles)}


# ----------------------------------------------------------------------------------------------- one run
def read_time(path) -> dict:
    """the /usr/bin/time record of a run (wall_s, peak_rss_kb, ...)"""
    out = {}
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                k, _, v = line.strip().partition("\t")
                if v:
                    try:
                        out[k] = float(v)
                    except ValueError:
                        pass
    return out


def check_rss(tm: dict, limit_gb: float = RSS_LIMIT_GB):
    gb = tm.get("peak_rss_kb", 0) / 1e6
    if gb > limit_gb:
        raise StopRun(f"peak RSS {gb:.1f} GB is above the stop limit of {limit_gb} GB")


def _session_members(sid: int) -> list[tuple[int, str, int]]:
    """(pid, state, resident kB) of every process of the session `sid`, read from /proc/<pid>/stat (state, session) and /proc/<pid>/statm (resident pages). No process
    is forked: the clean-up and the memory watchdog must work when `ps` cannot be started or stalls (the situations the watchdog exists for: a fork that fails, a
    swap storm). A process that ends while it is read is skipped."""
    out, page_kb = [], os.sysconf("SC_PAGE_SIZE") // 1024
    for ent in os.listdir("/proc"):
        if not ent.isdigit():
            continue
        try:
            with open(f"/proc/{ent}/stat") as fh:
                f = fh.read().rsplit(")", 1)[1].split()                  # after the command name: state, ppid, pgrp, session, ...
            if int(f[3]) != sid:
                continue
            with open(f"/proc/{ent}/statm") as fh:
                rss = int(fh.read().split()[1]) * page_kb
        except (OSError, ValueError, IndexError):
            continue
        out.append((int(ent), f[0], rss))
    return out


def group_rss_kb(sid: int) -> int:
    """summed resident set (kB) of every process of the session `sid` (the run's process group is its own session); a zombie holds no memory"""
    return sum(rss for _, state, rss in _session_members(sid) if not state.startswith("Z"))


def session_pids(sid: int) -> list[int]:
    """live (not zombie) processes of the session `sid`, which is the whole Liftoff tree: GNU time, liftoff, its pool workers, the cache wrapper and minimap2"""
    return [pid for pid, state, _ in _session_members(sid) if not state.startswith("Z")]


def kill_session(sid: int, *, term_grace_s: float = TERM_GRACE_S, kill_wait_s: float = KILL_WAIT_S, poll_s: float = 0.1) -> list[int]:
    """Empty the session `sid`: SIGTERM to its process group FIRST (one system call, before anything is listed: whatever happens to the listing, the group has been
    told), then to every member whatever its group, wait up to term_grace_s, SIGKILL the same way to what is left, wait up to kill_wait_s. Returns the pids that are
    STILL alive (an empty list is the normal result). For a session that is already empty the first signal goes to a group that no longer exists (ESRCH, ignored)."""
    if sid <= 1:
        raise ValueError(f"refusing to signal session {sid}")

    def hit(sig):
        try:
            os.killpg(sid, sig)
        except (ProcessLookupError, PermissionError):
            pass
        for pid in session_pids(sid):
            try:
                os.kill(pid, sig)
            except (ProcessLookupError, PermissionError):
                pass

    def wait(seconds):
        end = time.monotonic() + seconds
        while True:
            left = session_pids(sid)
            if not left or time.monotonic() >= end:
                return left
            time.sleep(poll_s)
    hit(signal.SIGTERM)
    left = wait(term_grace_s)
    if left:
        hit(signal.SIGKILL)
        left = wait(kill_wait_s)
    return left


_DEFER = {"on": False, "sig": None}          # set while the Liftoff tree is being created: a signal is remembered there, not raised (timed_watch)


@contextlib.contextmanager
def terminate_on_signals(signals=(signal.SIGTERM, signal.SIGHUP, signal.SIGINT)):
    """Inside the block a SIGTERM / SIGHUP / SIGINT raises SystemExit(128 + signal) in the main thread, so the `finally` blocks run (Python's default for SIGTERM
    and SIGHUP is to die at once, with the Liftoff session still running; the `timeout` of tools/rlock.sh sends SIGTERM). The old handlers are restored on exit.
    While _DEFER["on"] is set (the child is being created and its session is not known yet) the signal is only recorded in _DEFER["sig"]."""
    if threading.current_thread() is not threading.main_thread():
        yield
        return

    def handler(signum, frame):
        if _DEFER["on"]:
            _DEFER["sig"] = _DEFER["sig"] or signum
            return
        raise SystemExit(128 + signum)
    old = {s: signal.signal(s, handler) for s in signals}
    try:
        yield
    finally:
        for s, h in old.items():
            signal.signal(s, h)


@contextlib.contextmanager
def signals_held(signals=(signal.SIGTERM, signal.SIGHUP, signal.SIGINT)):
    """the clean-up of a session is not interrupted by a (second) signal: pending signals are delivered when the block ends"""
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    old = signal.pthread_sigmask(signal.SIG_BLOCK, signals)
    try:
        yield
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, old)


RC_TOO_LONG, RC_TOO_BIG = 124, 125


def timed_watch(cmd, time_file, log_file, *, cwd, hard_s: float, env: dict | None = None, rss_limit_gb: float = RSS_LIMIT_GB, poll_s: float = POLL_S,
                term_grace_s: float = TERM_GRACE_S, insession_margin_s: float = INSESSION_MARGIN_S, pass_fds: tuple = ()):
    """figures/_liftoff.py::_timed (own session under /usr/bin/time, killed after hard_s: returns 124) plus a watchdog: the session is killed while it runs when
    its summed RSS goes above rss_limit_gb (returns 125). The self-lift code checks the peak only after the run, which cannot stop a runaway one.
    Returns (return code, highest summed RSS seen in kB).

    The Liftoff tree must not outlive the driver (8 to 14 GB: two of them are the documented crash), so
      * the tree starts under `timeout -s KILL <hard_s + margin>` INSIDE its own session: even when the driver is SIGKILLed the session ends by itself, whatever its
        members do with SIGTERM (GNU timeout's `-k` would never fire once `time` has died of the first SIGTERM; the live driver sends SIGTERM and then SIGKILL at
        hard_s, earlier than this backstop);
      * a SIGTERM / SIGHUP / SIGINT to the driver raises SystemExit, and a `finally` kills the session whatever happened (a return, an exception, the signal); a signal
        that lands while the child is being created is remembered and raised as soon as the session is known (the handler must not block signals around Popen: the
        child would inherit the blocked mask);
      * the kill is SIGTERM, then SIGKILL of every member of the session, and the session is verified empty before this function returns or raises; it reads /proc and
        signals the group first, so it does not depend on `ps`;
      * `pass_fds` (the driver's run.lock) is inherited by `timeout`, `time`, the Liftoff main process and its pool workers, so the lock stays held while one of them
        lives; minimap2 and the cache wrapper do not inherit it (Liftoff starts them with close_fds), the stray guard covers them at the next call."""
    e = dict(os.environ)
    e["PATH"] = f"{L.LIFTOFF_ENV}/bin:" + e.get("PATH", "")
    e.setdefault("TMPDIR", "/mnt/linuxdisk/tmp")
    e.update(env or {})
    full = [TIMEOUT_BIN, "-s", "KILL", f"{hard_s + insession_margin_s:.3f}", TIME_BIN, "-f", L.TIME_FMT, "-o", str(time_file)] + [str(c) for c in cmd]
    peak, t0, p = 0, time.monotonic(), None
    with terminate_on_signals():
        try:
            with open(log_file, "w") as lf:
                _DEFER.update(on=True, sig=None)                    # a signal that lands while the child is being created is remembered, not raised: `p` is not assigned yet
                try:
                    p = subprocess.Popen(full, cwd=cwd, stdin=subprocess.DEVNULL, stdout=lf, stderr=subprocess.STDOUT, env=e, start_new_session=True,
                                         pass_fds=tuple(pass_fds))
                finally:
                    _DEFER["on"] = False
                if _DEFER["sig"]:
                    raise SystemExit(128 + _DEFER["sig"])           # now `p` is set: the `finally` below kills the session
                while True:
                    try:
                        rc = p.wait(timeout=max(0.05, min(poll_s, hard_s - (time.monotonic() - t0))))
                        # `timeout -s KILL` kills its whole group, itself included: a backstop that fired before the driver's own limit (a driver that was too slow to
                        # notice it) shows as a SIGKILLed leader. The run did outlive its limit: it is "too long", as `timeout` ending it would have been.
                        return (RC_TOO_LONG if rc == -signal.SIGKILL else rc), peak
                    except subprocess.TimeoutExpired:
                        pass
                    try:
                        peak = max(peak, group_rss_kb(p.pid))
                    except (OSError, ValueError, subprocess.SubprocessError):         # a poll that fails is skipped; it must not end the driver and leave the tree
                        pass
                    if peak / 1e6 > rss_limit_gb:
                        return RC_TOO_BIG, peak
                    if time.monotonic() - t0 >= hard_s:
                        return RC_TOO_LONG, peak
        finally:
            if p is not None:
                with signals_held():
                    left = kill_session(p.pid, term_grace_s=term_grace_s)
                    try:
                        p.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        pass
                    if left:
                        raise StopRun(f"the Liftoff session {p.pid} is not empty after SIGKILL (pids {left}): do not start another run until `ps --sid {p.pid}` is empty")


def find_stray_processes(e_root, proc: str = "/proc") -> list[tuple[int, str]]:
    """(pid, command line) of every process of THIS USER that is a liftoff, a minimap2 or the cache wrapper (one of the first two words of its command line) and
    that has `e_root` in its command line or as its working directory: what a killed driver may have left behind. Reads /proc; no pkill, no ps parsing."""
    roots = {str(e_root), os.path.realpath(e_root)}
    me, uid, hits = os.getpid(), os.getuid(), []
    for ent in os.listdir(proc):
        if not ent.isdigit() or int(ent) == me:
            continue
        try:
            if os.stat(f"{proc}/{ent}").st_uid != uid:
                continue
            with open(f"{proc}/{ent}/cmdline", "rb") as fh:
                argv = [a.decode(errors="replace") for a in fh.read().split(b"\0") if a]
            cwd = os.readlink(f"{proc}/{ent}/cwd")
        except OSError:
            continue
        if not argv or not ({os.path.basename(a) for a in argv[:2]} & STRAY_NAMES):
            continue
        text = " ".join(argv)
        if any(re.search(re.escape(r) + r"(?![\w.\-])", text) or cwd == r or cwd.startswith(r + os.sep) for r in roots):         # a path, not a sibling with the same prefix
            hits.append((int(ent), text[:300]))
    return sorted(hits)


STRAY_NAMES = {"liftoff", "minimap2", "mm2_cache.sh"}


def guard_no_stray(e_root):
    """a run starts only when no liftoff / minimap2 of this root is alive (an orphan of a killed driver; two of them at 8 to 14 GB are the documented crash)"""
    hits = find_stray_processes(e_root)
    if hits:
        raise StopRun(f"{len(hits)} liftoff / minimap2 process(es) of {e_root} are already running (pid {', '.join(str(p) for p, _ in hits[:6])}: "
                      f"{hits[0][1][:120]}): wait for them or stop them by pid, then call again")


class DriverLock:
    """a non-blocking exclusive flock on E_ROOT/run.lock for the life of one driver call: two drivers cannot overlap. `timeout`, `time`, the Liftoff main process and
    its pool workers inherit the file descriptor (timed_watch(pass_fds=(lock.fd,))), so the lock also stays held while one of them, left behind by a killed driver,
    lives. minimap2 and the cache wrapper do not inherit it (Liftoff starts them with close_fds): after a killed driver they are found by guard_no_stray, not by the lock."""

    def __init__(self, e_root):
        self.path = Path(e_root) / "run.lock"
        self.fd = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(self.path, os.O_RDWR | os.O_CREAT, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(fd)
            raise StopRun(f"{self.path} is held: another driver call, or the leftover Liftoff processes of a killed one (they inherit the lock), is running; "
                          f"look at `ps -eo pid,sid,etime,args | grep -E 'liftoff|minimap2'` and call again when they are gone") from None
        os.ftruncate(fd, 0)
        os.write(fd, f"{os.getpid()}\n".encode())
        self.fd = fd
        return self

    def __exit__(self, *exc):
        if self.fd is not None:
            os.close(self.fd)
            self.fd = None
        return False


def ensure_run_gffs(gff_path: str, assign: dict[str, str], out_dir) -> bool:
    """write the per-run reference GFFs (one pass over the annotation, a few seconds) unless they were written for this very assignment AND this very annotation
    file (size and mtime); True if written"""
    st = os.stat(gff_path)
    digest = hashlib.sha256(json.dumps([sorted(assign.items()), st.st_size, st.st_mtime_ns]).encode()).hexdigest()
    mark = Path(out_dir) / "assign.sha256"
    if mark.exists() and mark.read_text().strip() == digest and all((Path(out_dir) / lab / "in.gff").exists() for lab in set(assign.values())):
        return False
    write_run_gffs(gff_path, assign, str(out_dir))
    _atomic_write(mark, digest + "\n")
    return True


def plan_digest(plans: list[dict], hub_mode: str) -> str:
    """sha256 over the hub reading and every (run label, record ids) of the plans: what the runs of one E phase must all be cut from"""
    h = hashlib.sha256(hub_mode.encode())
    for run in sorted(plans, key=lambda r: r["label"]):
        h.update(f"|{run['label']}:{','.join(sorted(r['id'] for r in run['records']))}".encode())
    return h.hexdigest()


def _atomic_write(path, text: str):
    """a small state file is written whole or not at all (a kill in the middle leaves the old file, never half of the new one)"""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, p)


def freeze_plan(e_root, digest: str, write: bool = True):
    """the first call of the E phase records the plan digest; a later call with another plan (another hub reading, another pool table) is a stop.
    write=False only compares (the driver writes at the moment the first Liftoff process is about to start)."""
    p = Path(e_root) / "plan.digest"
    if p.exists():
        if p.read_text().strip() != digest:
            raise StopRun(f"{p} holds another plan digest: the runs of this E phase were cut from a different plan")
        return
    if write:
        _atomic_write(p, digest + "\n")


def check_plan_frozen(e_root, digest: str):
    """build: the plan must be the one the runs were cut from (nothing to compare when no run was made)"""
    p = Path(e_root) / "plan.digest"
    if p.exists() and p.read_text().strip() != digest:
        raise StopRun(f"the plan digest {digest[:12]} differs from {p.read_text().strip()[:12]} recorded by the first call of the E phase")


# Everything the run imports or executes: the two modules, the self-lift module and the figure library it imports, the lock wrapper of the heavy call, and the three test files
CODE_FILES = ("bench/hierarchy/expansion_table.py", "bench/hierarchy/liftoff_cross.py", "bench/hierarchy/test_expansion_table.py", "bench/hierarchy/test_supplement_e.py",
              "bench/hierarchy/test_supplement_e2.py", "figures/_liftoff.py", "figures/figlib.py", "tools/rlock.sh")


def resolve_commit(arg: str, repo=REPO) -> str:
    """the 40-hex commit that `arg` (an abbreviation or the hash) names in `repo`: `git rev-parse --verify <arg>^{commit}`. The run identity stores and compares THIS
    value, so an abbreviation in one call and the full hash (or another abbreviation) in the next are the same Amendment commit."""
    r = subprocess.run(["git", "-C", str(repo), "rev-parse", "--verify", f"{arg}^{{commit}}"], capture_output=True, text=True)
    full = r.stdout.strip()
    if r.returncode != 0 or not re.fullmatch(r"[0-9a-f]{40}", full):
        raise StopRun(f"{arg} is not a commit of {repo}")
    return full


CODE_COMMIT_PLACEHOLDER = "<CODE_COMMIT>"      # refused when still in the committed text; `<AMENDMENT_COMMIT>` in its command block is legitimate: a commit cannot contain its own hash
CODE_COMMIT_SENTENCE = re.compile(r"Commit\s+`([0-9a-f]{40})`\s+holds")


def check_committed_code(amendment_commit: str, repo=REPO, files=CODE_FILES, amendment_file: str | None = None, needles: dict | None = None) -> str:
    """the heavy run starts only from code that is committed: the Amendment commit is in the history of HEAD and every code file is tracked and equal to
    HEAD (compared by content, so `git update-index --assume-unchanged` hides nothing). With `amendment_file`, the text of that file AT the Amendment commit must
    have the heading AMENDMENT_MARK at the start of a line, not the literal `<CODE_COMMIT>` (`<AMENDMENT_COMMIT>` in the command block stays: a commit cannot contain
    its own hash), every hex digest of
    `needles` ({label: sha256}: the code files and the tools that would run, and on the first call the plan digest and the pinned inputs), and the sentence
    'Commit `<40-hex>` holds' naming a code commit that exists and holds every listed code file with the listed sha256, so the code that runs is the code the
    Amendment describes. Returns HEAD."""
    def git(*a):
        return subprocess.run(["git", "-C", str(repo), *a], capture_output=True, text=True)
    if git("cat-file", "-e", f"{amendment_commit}^{{commit}}").returncode != 0:
        raise StopRun(f"{amendment_commit} is not a commit of {repo}")
    if git("merge-base", "--is-ancestor", amendment_commit, "HEAD").returncode != 0:
        raise StopRun(f"HEAD does not contain the Amendment commit {amendment_commit}")
    for f in files:
        if git("ls-files", "--error-unmatch", f).returncode != 0:
            raise StopRun(f"{f} is not tracked: commit the code before the run")
        if git("diff", "--quiet", "HEAD", "--", f).returncode != 0 or git("hash-object", "--", f).stdout.strip() != git("rev-parse", f"HEAD:{f}").stdout.strip():
            raise StopRun(f"{f} differs from HEAD: commit the code before the run")
    if amendment_file is not None:
        shown = git("show", f"{amendment_commit}:{amendment_file}")
        if shown.returncode != 0:
            raise StopRun(f"{amendment_file} is not in the Amendment commit {amendment_commit}")
        text = shown.stdout
        if not re.search(r"^" + re.escape(AMENDMENT_MARK) + r"\b", text, re.M):
            raise StopRun(f"{amendment_file} at {amendment_commit} has no '{AMENDMENT_MARK}' heading at the start of a line: that commit does not hold the Amendment")
        if CODE_COMMIT_PLACEHOLDER in text:
            raise StopRun(f"the Amendment at {amendment_commit} still contains the placeholder {CODE_COMMIT_PLACEHOLDER}: name the code commit and commit again before the run")
        for label, hexdigest in (needles or {}).items():
            if hexdigest not in text:
                raise StopRun(f"the Amendment at {amendment_commit} does not list the sha256 of {label} ({hexdigest[:16]}...): what would run or be frozen is not what it describes")
        m = CODE_COMMIT_SENTENCE.search(text)
        if not m:
            raise StopRun(f"the Amendment at {amendment_commit} does not name its code commit as 'Commit `<40-hex>` holds' (section 1): the code it lists cannot be traced")
        code_commit = m.group(1)
        if git("cat-file", "-e", f"{code_commit}^{{commit}}").returncode != 0:
            raise StopRun(f"the code commit {code_commit} named by the Amendment is not a commit of {repo}")
        for f in files:
            want = (needles or {}).get(f)
            if want is None:
                continue
            blob = subprocess.run(["git", "-C", str(repo), "show", f"{code_commit}:{f}"], capture_output=True)
            if blob.returncode != 0:
                raise StopRun(f"the code commit {code_commit[:12]} named by the Amendment does not contain {f}")
            if hashlib.sha256(blob.stdout).hexdigest() != want:
                raise StopRun(f"the code commit {code_commit[:12]} named by the Amendment holds {f} with another sha256 than the one the Amendment lists")
    return git("rev-parse", "HEAD").stdout.strip()


def code_identity(repo=REPO, files=CODE_FILES) -> dict:
    return {f: X._sha(os.path.join(str(repo), f)) for f in files}


def _sha_file(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def tool_identity(liftoff_bin, minimap2, time_bin=TIME_BIN) -> dict:
    """sha256 of the executables a run uses: the Liftoff entry script, the minimap2 it calls (the real one behind the cache wrapper) and GNU time"""
    return {name: (_sha_file(path) if os.path.isfile(path) else "absent") for name, path in (("liftoff", liftoff_bin), ("minimap2", minimap2), ("gnu_time", time_bin))}


def command_fingerprint(liftoff_bin, reference, minimap2) -> dict:
    """what the Liftoff command line is made of (the target is a placeholder: it differs by outgroup)"""
    return {"liftoff_bin": str(liftoff_bin), "reference": str(reference), "minimap2": str(minimap2), "cmd": liftoff_command("TARGET", reference, minimap2, liftoff_bin),
            "params": L.params_key()}


def input_paths(gorilla_fasta, gorilla_gff, lift_root, outgroups) -> dict:
    """name -> path of every input of Amendment 2 section 3 (17 files): the gorilla assembly and its .fai, its annotation, its self-lift table and record table,
    and for every outgroup its assembly, .fai, k15 w10 index and own self-lift table"""
    d = {"gorilla_fasta": str(gorilla_fasta), "gorilla_fai": str(gorilla_fasta) + ".fai", "gorilla_gff": str(gorilla_gff),
         "gorilla_selflift": str(lift_root) + "/gorilla/liftoff_loci.tsv", "gorilla_records": str(lift_root) + "/gorilla/records.tsv"}
    for k, t in sorted(outgroups.items()):
        d.update({f"{k}_fasta": t["genome"], f"{k}_fai": t["genome"] + ".fai", f"{k}_mmi": t["mmi"], f"{k}_selflift": t["selflift"]})
    return d


def pin_inputs(paths: dict, deep: bool = True) -> dict:
    """{name: {path, size, mtime_ns, sha256}}: a file below INPUT_BIG_BYTES carries its sha256 (when `deep`, else only below SHALLOW_SHA_BYTES), a bigger one only
    its size and mtime (its sha256 is held by Gate 0)"""
    out = {}
    for name in sorted(paths):
        path = str(paths[name])
        try:
            st = os.stat(path)
        except OSError as e:
            raise StopRun(f"input {name}: {path} cannot be read ({e.strerror})") from e
        rec = {"path": path, "size": st.st_size, "mtime_ns": st.st_mtime_ns, "sha256": None}
        if st.st_size < INPUT_BIG_BYTES and (deep or st.st_size < SHALLOW_SHA_BYTES):
            rec["sha256"] = _sha_file(path)
        out[name] = rec
    return out


def compare_pins(recorded: dict, current: dict) -> list[str]:
    """what differs between the pins of the first call and the pins of this one (empty = nothing); sha256 is compared when both sides have it"""
    problems = []
    for name in sorted(set(recorded) | set(current)):
        r, c = recorded.get(name), current.get(name)
        if r is None:
            problems.append(f"input {name} ({c['path']}) was not pinned by the first call")
            continue
        if c is None:
            problems.append(f"input {name} ({r['path']}) is no longer an input of this call")
            continue
        what = []
        if r["path"] != c["path"]:
            what.append(f"path {r['path']} -> {c['path']}")
        if r["size"] != c["size"]:
            what.append(f"size {r['size']} -> {c['size']}")
        if r["mtime_ns"] != c["mtime_ns"]:
            what.append(f"mtime {r['mtime_ns']} -> {c['mtime_ns']}")
        if r.get("sha256") and c.get("sha256") and r["sha256"] != c["sha256"]:
            what.append("sha256")
        if what:
            problems.append(f"input {name} ({c['path']}) changed since the first call: {', '.join(what)}")
    return problems


IDENTITY_NAMES = {"code_sha256": "code", "amendment_commit": "Amendment commit", "tools": "tools", "command": "command line"}


def load_identity(e_root) -> dict | None:
    p = Path(e_root) / "run_identity.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except ValueError as e:
        raise StopRun(f"{p} is not valid JSON ({e}): the identity of the E phase is lost; do not guess, report") from e


def check_identity(path, have: dict, want: dict, pins: dict | None):
    """StopRun unless every key of `want` equals the recorded one and (when `pins` is given) no pinned input changed"""
    bad = [IDENTITY_NAMES.get(k, k) for k in want if have.get(k) != want[k]]
    if bad:
        raise StopRun(f"the run identity {path} records other {' and other '.join(bad)} than this call: the runs of one E phase use one code, one set of tools and one command")
    if pins is not None:
        problems = compare_pins(have.get("inputs", {}), pins)
        if problems:
            raise StopRun(f"the run identity {path}: " + "; ".join(problems[:4]) + (f" (and {len(problems) - 4} more)" if len(problems) > 4 else ""))


def freeze_identity(e_root, identity: dict, amendment_commit: str, write: bool = True, extra: dict | None = None, pins: dict | None = None):
    """the first call records the sha256 of the code files, the Amendment commit (and `extra`: the tools and the command, and the pinned inputs `pins`); a later
    call with other code, another Amendment, other tools, another command or a changed input is a stop. write=False only compares (the driver writes at the moment
    the first Liftoff process is about to start)."""
    p = Path(e_root) / "run_identity.json"
    want = {"amendment_commit": amendment_commit, "code_sha256": identity, **(extra or {})}
    have = load_identity(e_root)
    if have is not None:
        check_identity(p, have, want, pins)
        return
    if write:
        rec = dict(want)
        if pins is not None:
            rec["inputs"] = pins
        _atomic_write(p, json.dumps(rec, indent=1, sort_keys=True) + "\n")


def verify_identity(e_root, *, code: dict, tools: dict, command: dict, pins: dict):
    """build: run_identity.json must exist and describe the code, the tools, the command and the inputs of this call"""
    p = Path(e_root) / "run_identity.json"
    have = load_identity(e_root)
    if have is None:
        raise StopRun(f"{p} is missing: nothing in {e_root} was made by a `run` call of the Amendment (or its identity was deleted)")
    check_identity(p, have, {"code_sha256": code, "tools": tools, "command": command}, pins)


def _key(cmd: list[str], gff: Path) -> dict:
    h = hashlib.sha256(gff.read_bytes()).hexdigest()[:16] if gff.exists() else "absent"
    return {"cmd": cmd, "in_gff_sha256": h, "params": L.params_key()}


def parent_ids(gff) -> set[str]:
    """ids of the parent records (gene / pseudogene without a Parent) of a GFF3: what Liftoff lifts"""
    out = set()
    with open(gff) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in L.PARENT_TYPES or _PARENT.search(f[8]):
                continue
            m = _ID.search(f[8])
            if m:
                out.add(m.group(1))
    return out


def check_accounted(d: Path, gff: Path):
    """Liftoff does not look at the exit status of the minimap2 it runs: a minimap2 that died or was killed leaves a short SAM and Liftoff still exits 0 with a
    short lifted.gff3. In every real self-lift (179,088 records in four genomes) each parent record of the input is either lifted (a primary or a copy row of
    loci.tsv) or listed in unmapped.txt, so a record that is neither is the signature of a failed alignment: the run is not a result (StopRun, no done.key)."""
    up = d / "unmapped.txt"
    if not up.exists():
        raise StopRun(f"{up} is missing: Liftoff writes one for every run, so this run did not end the way Liftoff ends (see {d / 'liftoff.log'})")
    parents = parent_ids(gff)
    unmapped = {x.strip() for x in up.read_text().splitlines() if x.strip()}
    with open(d / "loci.tsv") as fh:
        lifted = {X.source_of(r["feature_id"], parents) for r in csv.DictReader(fh, delimiter="\t")}
    missing = sorted(parents - lifted - unmapped)
    if missing:
        (d / "incomplete").write_text(json.dumps({"missing": missing, "parents": len(parents)}) + "\n")
        raise StopRun(f"liftoff exited 0 but {len(missing)} of {len(parents)} parent record(s) of {gff} are neither lifted nor listed in unmapped.txt "
                      f"({', '.join(missing[:6])}{' ...' if len(missing) > 6 else ''}): the minimap2 under Liftoff probably died, and Liftoff does not look at its exit status; "
                      f"see {d / 'liftoff.log'} (lifted.gff3 and loci.tsv are kept, there is no done.key)")


def run_shard(run_dir, record_ids, *, liftoff_bin: str, target: str, reference: str, minimap2: str, hard_s: int = HARD_LIMIT_S, env: dict | None = None,
              rss_limit_gb: float = RSS_LIMIT_GB, poll_s: float = POLL_S, term_grace_s: float = TERM_GRACE_S, insession_margin_s: float = INSESSION_MARGIN_S,
              pass_fds: tuple = (), timing: dict | None = None) -> dict:
    """Run Liftoff for the records of one shard, in `run_dir` (which holds in.gff). Returns {"status": "done" | "skipped" | "toolong", ...}. A run that is
    killed at the hard limit leaves a `toolong` file and no done.key (the driver splits it); a run whose memory goes above the limit is killed while it
    runs; any other failure raises StopRun (no retry). A run is `done` only if Liftoff exited 0 AND every parent record of in.gff is lifted or listed in
    unmapped.txt (check_accounted). The files of an earlier attempt are removed before the launch; `timing["elapsed_s"]` is set whatever happens."""
    d = Path(run_dir)
    gff = d / "in.gff"
    if not gff.exists():
        raise StopRun(f"{gff} is missing: write the run GFFs first")
    cmd = liftoff_command(target, reference, minimap2, liftoff_bin)
    key = _key(cmd, gff)
    kpath = d / "done.key"
    if _is_done(d, key):
        return {"status": "skipped"}
    for p in (kpath, d / "toolong", d / "rss_stop", d / "incomplete", d / "loci.tsv", d / "time.txt", d / "lifted.gff3", d / "lifted.gff3.gz", d / "unmapped.txt",
              Path(str(gff) + "_db")):
        p.unlink(missing_ok=True)
    shutil.rmtree(d / "inter", ignore_errors=True)
    (d / "types.txt").write_text(types_text())
    t_launch = time.monotonic()
    try:
        try:
            rc, peak_kb = timed_watch(cmd, d / "time.txt", d / "liftoff.log", cwd=d, hard_s=hard_s, env=env, rss_limit_gb=rss_limit_gb, poll_s=poll_s,
                                      term_grace_s=term_grace_s, insession_margin_s=insession_margin_s, pass_fds=pass_fds)
        finally:
            elapsed = time.monotonic() - t_launch
            if timing is not None:
                timing["elapsed_s"] = elapsed
        tm = read_time(d / "time.txt")
        if rc == RC_TOO_BIG:
            (d / "rss_stop").write_text(f"killed at a summed RSS of {peak_kb / 1e6:.1f} GB (limit {rss_limit_gb} GB), {len(record_ids)} records\n")
            raise StopRun(f"RSS {peak_kb / 1e6:.1f} GB is above the stop limit of {rss_limit_gb} GB in {d}: the run was killed")
        if rc == RC_TOO_LONG:
            (d / "toolong").write_text(f"killed after {hard_s} s at {len(record_ids)} records\n")
            return {"status": "toolong", "records": len(record_ids), "wall_s": tm.get("wall_s"), "elapsed_s": elapsed, "peak_summed_kb": peak_kb}
        if rc != 0 or not (d / "lifted.gff3").exists():
            raise StopRun(f"liftoff failed (exit {rc}) in {d}: see {d / 'liftoff.log'}")
        check_rss(tm, rss_limit_gb)
        n = L.parse_lifted(d / "lifted.gff3", d / "loci.tsv")
        check_accounted(d, gff)
        with open(d / "lifted.gff3", "rb") as fi, gzip.open(d / "lifted.gff3.gz", "wb", compresslevel=4) as fo:
            shutil.copyfileobj(fi, fo)
        (d / "lifted.gff3").unlink()
        Path(str(gff) + "_db").unlink(missing_ok=True)
        _atomic_write(kpath, json.dumps(key))
        return {"status": "done", "features": n, "wall_s": tm.get("wall_s"), "elapsed_s": elapsed, "peak_rss_kb": tm.get("peak_rss_kb"), "peak_summed_kb": peak_kb}
    finally:
        shutil.rmtree(d / "inter", ignore_errors=True)


# ----------------------------------------------------------------------------------------------- the schedule
GROUPS = ("pos", "neg", "verdict")


def schedule_units(verdict_plan: list[dict], pos_plan: list[dict], neg_plan: list[dict], outgroups=("H", "C", "O")) -> list[tuple[str, str]]:
    """(outgroup, run label) in execution order: the positive-control runs (the LRPAP1 and fragment components), then the planted-negative runs, then the
    ranks of the verdict set; within a group outgroup by outgroup. Both controls are lifted and checked before any verdict-set record is lifted."""
    units = []
    for plan in (pos_plan, neg_plan, verdict_plan):
        for o in outgroups:
            units += [(o, r["label"]) for r in plan]
    return units


def group_map(verdict_plan: list[dict], pos_plan: list[dict], neg_plan: list[dict]) -> dict[str, str]:
    """run label -> group (pos, neg or verdict)"""
    out = {}
    for g, plan in (("pos", pos_plan), ("neg", neg_plan), ("verdict", verdict_plan)):
        for r in plan:
            out[r["label"]] = g
    return out


def predict_wall(n_records: int, per_record_s: float = PER_RECORD_S, startup_s: float = STARTUP_S, margin: float = 1.25) -> float:
    return startup_s + n_records * per_record_s * margin


def plan_calls(units: list, predicted: dict, budget_s: float = CALL_BUDGET_S, group_of: dict | None = None) -> list[list]:
    """cut the unit list into bounded calls: a unit starts only if its predicted time fits what is left of the call; the first unit of a call always runs;
    with `group_of` a call ends where the group changes (the controls have calls of their own)"""
    calls, cur, used, cur_g = [], [], 0.0, None
    for u in units:
        t = predicted[u]
        g = group_of[u[1]] if group_of else None
        if cur and (used + t > budget_s or g != cur_g):
            calls.append(cur)
            cur, used = [], 0.0
        cur.append(u)
        used += t
        cur_g = g
    if cur:
        calls.append(cur)
    return calls


# ----------------------------------------------------------------------------------------------- the budgeted driver
def load_splits(e_root, species) -> list[str]:
    p = Path(e_root) / "lift" / species / "splits.json"
    if not p.exists():
        return []
    try:
        return json.loads(p.read_text())
    except ValueError as e:
        raise StopRun(f"{p} is not valid JSON ({e}): the shards that outlived the hard limit are lost; do not guess, report") from e


def save_splits(e_root, species, splits: list[str]):
    _atomic_write(Path(e_root) / "lift" / species / "splits.json", json.dumps(sorted(set(splits))))


def effective_shards(shards, by_id, splits) -> list[tuple[str, list[str]]]:
    """the planned shards of a run with every shard that outlived the hard limit replaced by its two halves (again if a half outlives it)"""
    todo, out = list(shards), []
    while todo:
        lab, ids = todo.pop(0)
        if lab in splits:
            try:
                todo[0:0] = X.split_shard(lab, ids, by_id)
            except X.NoSplit as e:
                raise StopRun(str(e)) from e
        else:
            out.append((lab, ids))
    return out


def _install_gff(e_root, run_label, shard_label, ids, d: Path):
    """in.gff of a shard: the run's reference GFF, or its records of this shard when the run was blocked or split"""
    src = Path(e_root) / "gff" / run_label / "in.gff"
    if not src.exists():
        raise StopRun(f"{src} is missing: write the run GFFs first (`gffs`)")
    if shard_label == run_label:
        data = src.read_bytes()
        if not (d / "in.gff").exists() or (d / "in.gff").read_bytes() != data:
            (d / "in.gff").write_bytes(data)
        return
    write_run_gffs(str(src), {rid: shard_label for rid in ids}, str(d.parent))


def shard_source_gff(e_root, run_label, shard_label, d: Path) -> Path:
    """the file a shard's in.gff is (or is cut from): the run's reference GFF; a split or blocked shard is cut from it, so its own in.gff is what there is"""
    return Path(e_root) / "gff" / run_label / "in.gff" if shard_label == run_label else Path(d) / "in.gff"


def _is_done(d: Path, key: dict | None = None) -> bool:
    """a run is finished iff done.key and loci.tsv exist and, when `key` is given, done.key holds exactly that key (the command with its Liftoff binary, target,
    reference and minimap2, the in.gff, the parameters): output made by another command or from another in.gff is not a finished run"""
    kp = Path(d) / "done.key"
    if not (kp.exists() and (Path(d) / "loci.tsv").exists()):
        return False
    if key is None:
        return True
    try:
        return json.loads(kp.read_text()) == key
    except ValueError:
        return False


class RunSpec:
    """what a finished run of this root must have been made with (used by the gate and by build to refuse foreign output)"""

    def __init__(self, e_root, outgroups, reference, liftoff_bin, minimap2):
        self.e_root, self.outgroups, self.reference, self.liftoff_bin, self.minimap2 = e_root, outgroups, reference, liftoff_bin, minimap2

    def key(self, outgroup: str, run_label: str, shard_label: str) -> dict:
        t = self.outgroups[outgroup]
        target = str(Path(self.e_root) / "lift" / t["species"] / "genome.fa")
        d = run_dir_of(outgroup, shard_label, self.e_root, self.outgroups)
        return _key(liftoff_command(target, self.reference, self.minimap2, self.liftoff_bin), shard_source_gff(self.e_root, run_label, shard_label, d))


class WallLedger:
    """E_ROOT/wall.json: the wall seconds of every Liftoff process the driver started, per outgroup, over all the calls of the phase (stop rule 6)"""

    def __init__(self, e_root):
        self.path = Path(e_root) / "wall.json"
        self.data: dict = {}
        if self.path.exists():
            try:
                self.data = json.loads(self.path.read_text())
            except ValueError as e:
                raise StopRun(f"{self.path} is not valid JSON ({e}): the Liftoff wall of the phase is lost; do not guess, report") from e

    def entries(self, outgroup: str) -> list[dict]:
        return list(self.data.get(outgroup, []))

    def total(self, outgroup: str) -> float:
        return sum(e["elapsed_s"] for e in self.entries(outgroup))

    def add(self, outgroup: str, label: str, elapsed_s: float, status: str):
        self.data.setdefault(outgroup, []).append(dict(label=label, elapsed_s=round(float(elapsed_s), 3), status=status))
        _atomic_write(self.path, json.dumps(self.data, indent=1, sort_keys=True) + "\n")


def predicted_wall_by_outgroup(units, plan_by_label, per_record_s: float = PER_RECORD_S, startup_s: float = STARTUP_S) -> dict:
    """the predicted Liftoff wall of every outgroup: the sum of predict_wall over all its units of the plan (the registered 1,669 s each)"""
    out: dict = collections.defaultdict(float)
    for o, label in units:
        out[o] += predict_wall(len(plan_by_label[label]["records"]), per_record_s, startup_s)
    return dict(out)


def run_budgeted(units, plan_by_label, *, e_root, outgroups, reference, liftoff_bin, minimap2, budget_s=CALL_BUDGET_S, hard_s=HARD_LIMIT_S, env=None,
                 per_record_s=PER_RECORD_S, startup_s=STARTUP_S, rss_limit_gb=RSS_LIMIT_GB, ensure_links=True, clock=time.monotonic,
                 group_of=None, gate=None, call_limit_s=CALL_LIMIT_S, log=None, t0=None, stray_guard=True, on_first_launch=None, pass_fds=(),
                 term_grace_s=TERM_GRACE_S, insession_margin_s=INSESSION_MARGIN_S, poll_s=POLL_S, wall_pred_s=None, wall_factor=None) -> dict:
    """One bounded call. Units run in order; a unit starts only if its predicted time fits what is left of the call (the first unit always runs). A shard
    killed at the hard limit is recorded in splits.json (the next call runs its halves; no record-free cut = StopRun); any failure raises StopRun.
    With `group_of` (label -> group) a call serves ONE group: the group of the first unit that still has work, `gate(group)` is called before that unit
    starts (it raises StopRun to refuse), and the call ends where the group changes. Units that are already done (by THIS command, from THIS in.gff: see
    _is_done) are skipped without a gate.
    No run outlives `call_limit_s` seconds from the start of its call. A run killed at its hard limit is too long and its shard is split; a run that is not the
    first of its call and is killed at the end of the call instead (it had less than the hard limit) is not too long, the next call runs it again ("deadline").
    Returns {"done", "pending", "split", "deadline", "group", "exit_code", "wall_s"}: exit_code 75 = run the same command again.

    `t0` is the start of the call on `clock` (the driver passes the start of its process, so the start-up work is charged to the call); `stray_guard` refuses to
    start a run while a liftoff / minimap2 of this root is alive; `on_first_launch()` is called once, immediately before the first Liftoff process of the call
    is started (the driver freezes the plan and the identity there, so a call that cannot start leaves nothing frozen); `pass_fds` are inherited by the tree.
    `wall_pred_s` ({outgroup: predicted wall in s}) switches on stop rule 6: the wall of every Liftoff process is added to E_ROOT/wall.json, and no run starts for
    an outgroup whose sum is above `wall_factor` (default WALL_FACTOR) times its prediction (StopRun)."""
    t0 = clock() if t0 is None else t0
    done, ran_any, group = [], False, None
    announced = False
    log = log or (lambda msg: None)
    factor = WALL_FACTOR if wall_factor is None else wall_factor
    ledger = WallLedger(e_root) if wall_pred_s is not None else None

    def finish(i, exit_code, split=None, deadline=None, pending=None):
        return {"done": done, "pending": list(units[i:]) if pending is None else pending, "split": split, "deadline": deadline, "group": group, "exit_code": exit_code,
                "wall_s": {o: round(ledger.total(o)) for o in sorted(wall_pred_s)} if ledger is not None else {}}

    def check_wall(o):
        lim = (wall_pred_s or {}).get(o)
        if ledger is None or lim is None:
            return
        total, limit = ledger.total(o), factor * lim
        if total > limit:
            raise StopRun(f"time budget of {outgroups[o]['species']}: its Liftoff runs took {total:.0f} s in all, more than {factor:g} x the predicted {lim:.0f} s "
                          f"({limit:.0f} s); stop and report to the user (what finished stays valid)")
    for i, (o, label) in enumerate(units):
        run = plan_by_label[label]
        t = outgroups[o]
        if ensure_links:
            target = str(ensure_target_links(o, e_root, outgroups) / "genome.fa")
        else:
            target = t["genome"]
        cmd = liftoff_command(target, reference, minimap2, liftoff_bin)
        by_id = {r["id"]: r for r in run["records"]}
        shards = effective_shards(run["shards"], by_id, load_splits(e_root, t["species"]))
        todo = []
        for sl, ids in shards:
            d = run_dir_of(o, sl, e_root, outgroups)
            if _is_done(d, _key(cmd, shard_source_gff(e_root, label, sl, d))):
                continue
            if _is_done(d):
                log(f"{o} {sl}: done.key is not the key of this call's command and in.gff: the run is made again")
            todo.append((sl, ids))
        if not todo:
            done.append((o, label))
            continue
        if group_of:
            if group is None:
                group = group_of[label]
                if gate is not None:
                    gate(group)
                    log(f"gate({group}) passed")
            elif group_of[label] != group:
                log(f"call ends: {o} {label} belongs to another group than {group}")
                return finish(i, 75)
        pred = sum(predict_wall(len(ids), per_record_s, startup_s) for _, ids in todo)
        if ran_any and clock() - t0 + pred > budget_s:
            return finish(i, 75)
        for sl, ids in todo:
            check_wall(o)
            if stray_guard:
                guard_no_stray(e_root)
            if not announced:
                announced = True
                if on_first_launch is not None:
                    on_first_launch()
            d = run_dir_of(o, sl, e_root, outgroups)
            d.mkdir(parents=True, exist_ok=True)
            _install_gff(e_root, label, sl, ids, d)
            first_run = not ran_any
            hard_eff = min(hard_s, max(1.0, call_limit_s - (clock() - t0)))
            timing, t_launch = {}, time.monotonic()
            try:
                res = run_shard(d, ids, liftoff_bin=liftoff_bin, target=target, reference=reference, minimap2=minimap2, hard_s=hard_eff, env=env,
                                rss_limit_gb=rss_limit_gb, poll_s=poll_s, term_grace_s=term_grace_s, insession_margin_s=insession_margin_s, pass_fds=pass_fds,
                                timing=timing)
            except BaseException as e:                                  # a stop, a signal: the time the process took still counts
                if ledger is not None and timing.get("elapsed_s", 0) > 0:
                    ledger.add(o, sl, timing["elapsed_s"], "stopped" if isinstance(e, StopRun) else "interrupted")
                raise
            ran_any = True
            if ledger is not None and res["status"] != "skipped":
                ledger.add(o, sl, res.get("elapsed_s", timing.get("elapsed_s", 0.0)), res["status"])
            log(f"{o} {sl} [{group or '-'}]: {res['status']}, {len(ids)} records, wall {res.get('wall_s')} s, elapsed {res.get('elapsed_s', 0.0):.0f} s, "
                f"peak RSS {(res.get('peak_rss_kb') or 0) / 1e6:.1f} GB, summed {(res.get('peak_summed_kb') or 0) / 1e6:.1f} GB")
            if res["status"] == "toolong" and hard_eff < hard_s and not first_run:
                (d / "toolong").unlink(missing_ok=True)
                check_wall(o)
                return finish(i, 75, deadline=(o, label))
            if res["status"] == "toolong":
                by = {r["id"]: r for r in run["records"]}
                try:
                    X.split_shard(sl, ids, by)
                except X.NoSplit as e:
                    raise StopRun(f"{o} {sl}: killed at {hard_s} s and {e}") from e
                save_splits(e_root, t["species"], load_splits(e_root, t["species"]) + [sl])
                check_wall(o)
                return finish(i, 75, split=(o, label))
            check_wall(o)                                               # stop rule 6 right after the run that passes the limit (the state on disk is complete)
        done.append((o, label))
    return finish(len(units), 0, pending=[])


# ----------------------------------------------------------------------------------------------- work dirs (created only by the authorised run)
def _usable(path) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def check_liftoff_binary(path):
    """the Liftoff entry script exists and can be executed"""
    if not (os.path.isfile(path) and os.access(path, os.X_OK)):
        raise StopRun(f"the Liftoff binary {path} does not exist or is not executable")


def check_supervisors():
    """the two programs the Liftoff tree is started under, `timeout` (the in-session backstop) and GNU `time` (peak RSS and wall), exist"""
    for b in (TIMEOUT_BIN, TIME_BIN):
        if not (os.path.isfile(b) and os.access(b, os.X_OK)):
            raise StopRun(f"{b} does not exist or is not executable: the Liftoff tree is started under `timeout -s KILL ... time ...`")


def make_lift_dir(e_root) -> Path:
    """E_ROOT/lift, where the minimap2 cache wrapper and the per-outgroup work directories live (created before the wrapper is written into it)"""
    d = Path(e_root) / "lift"
    d.mkdir(parents=True, exist_ok=True)
    return d


def check_reference(fasta):
    """the gorilla assembly Liftoff reads: it and its .fai exist, and the .fai is not older (pyfaidx rewrites an index older than its FASTA)"""
    fai = str(fasta) + ".fai"
    for f in (str(fasta), fai):
        if not _usable(f):
            raise StopRun(f"{f} is missing or empty")
    if os.stat(fai).st_mtime < os.stat(fasta).st_mtime:
        raise StopRun(f"{fai} is older than {fasta}: pyfaidx would rewrite it")


def check_outgroup_files(t: dict):
    """the three files behind the links of an outgroup exist, are not empty, and nothing will rewrite them: the .fai is not older than the FASTA (pyfaidx would
    rebuild it THROUGH the link) and the .mmi is the k15 w10 index of this very FASTA (Liftoff builds `<target>.mmi` through a dangling link, 7.7 GB, when it is absent)"""
    fasta, fai, mmi = t["genome"], t["genome"] + ".fai", t["mmi"]
    for f in (fasta, fai, mmi):
        if not _usable(f):
            raise StopRun(f"{f} is missing or empty")
    if os.stat(fai).st_mtime < os.stat(fasta).st_mtime:
        raise StopRun(f"{fai} is older than {fasta}: pyfaidx would rewrite it through the link")
    try:
        w, k, b, n, flag = L._mmi_header(mmi)
    except (ValueError, OSError, struct.error) as e:
        raise StopRun(f"{mmi} is not a minimap2 index ({e})") from e
    with open(fai) as fh:
        n_seq = sum(1 for _ in fh)
    if (w, k, flag) != (10, 15, 0) or n != n_seq:
        raise StopRun(f"{mmi} is not the Liftoff index of {fasta}: the index has w{w} k{k} flag {flag} and {n} sequence(s), the registered one is w10 k15 flag 0 with {n_seq}")


def ensure_target_links(key: str, e_root, outgroups: dict | None = None) -> Path:
    """E_ROOT/lift/<species>/genome.fa, .fai and .mmi as links to the outgroup FASTA, its index and its Liftoff index: Liftoff reuses `<target>.mmi`
    when it exists (align_features.build_minimap2_index). The links point at files that are read, never written; the files are checked first
    (check_outgroup_files) and an existing link must be exactly the registered one."""
    t = (outgroups or OUTGROUPS)[key]
    check_outgroup_files(t)
    d = Path(e_root) / "lift" / t["species"]
    srcs = ((t["genome"], "genome.fa"), (t["genome"] + ".fai", "genome.fa.fai"), (t["mmi"], "genome.fa.mmi"))
    for src, name in srcs:
        dst = d / name
        if os.path.lexists(dst) and (not dst.is_symlink() or os.path.realpath(dst) != os.path.realpath(src)):
            raise StopRun(f"{dst} is not a link to {src}: it is stale; remove that link by hand if you are sure")
    d.mkdir(parents=True, exist_ok=True)
    for src, name in srcs:
        dst = d / name
        if not os.path.lexists(dst):
            dst.symlink_to(src)
    return d


def run_dir_of(key: str, label: str, e_root, outgroups: dict | None = None) -> Path:
    return Path(e_root) / "lift" / (outgroups or OUTGROUPS)[key]["species"] / "runs" / label


def plan_from_real_inputs(hub_mode: str = "nodes", lift_root: str | None = None, spec: dict | None = None, n_neg: int = N_PLANTED) -> dict:
    """The plans of the gorilla self-lift (input-only). called = the non-hub components with at least two such records (the registered called set);
    control = the LRPAP1 component and the fragment component, verdict = the rest of the called set; plan = the rank plan of the verdict set, pos = that of
    the control components (labels pos1..), ctl = the two runs of the planted negative control, called_plan = the rank plan of the whole called set (the
    numbers quoted in section 7). The ids of the records to be lifted are checked for the collision that would misread a copy row (X.check_no_id_collisions).
    `lift_root` None = X.LIFT_ROOT as it is when the function is called."""
    lift_root = X.LIFT_ROOT if lift_root is None else lift_root
    rows = X.read_selflift(lift_root + "/gorilla/liftoff_loci.tsv")
    recs = X.read_records(lift_root + "/gorilla/records.tsv")
    pool = X.build_pool(rows, recs)
    called = X.called_components(pool, hub_mode)
    ctrl = X.lrpap_control(pool, rows, spec)
    verdict, control = X.split_control(called, ctrl)
    neg = X.plant_negative(rows, recs, pool, n_records=n_neg)
    w = dict(hub_mode=hub_mode, rows=rows, recs=recs, pool=pool, called=called, ctrl=ctrl, verdict=verdict, control=control, neg=neg,
             called_plan=X.rank_plan(called), plan=X.rank_plan(verdict), pos=X.rank_plan(control, label_prefix="pos"), ctl=X.rank_plan(neg, label_prefix="ctl"))
    X.check_no_id_collisions([r["id"] for k in ("plan", "pos", "ctl") for run in w[k] for r in run["records"]])
    return w


# ----------------------------------------------------------------------------------------------- Gate 7 from the finished runs
def shard_pairs(plan: list[dict], e_root, outgroups, key: str) -> list[tuple[str, str]]:
    """(run label, shard label) of the shards that make up `plan` in outgroup `key`: the planned shards, each replaced by its halves where splits.json says so"""
    splits = load_splits(e_root, outgroups[key]["species"])
    out = []
    for run in plan:
        by_id = {r["id"]: r for r in run["records"]}
        out += [(run.get("label"), lab) for lab, _ in effective_shards(run["shards"], by_id, splits)]
    return out


def shard_labels(plan: list[dict], e_root, outgroups, key: str) -> list[str]:
    return [lab for _, lab in shard_pairs(plan, e_root, outgroups, key)]


def finished_shard_labels(plan: list[dict], e_root, outgroups, key: str, runspec: RunSpec | None = None) -> list[str]:
    """the shard labels of `plan` in outgroup `key`; every one must be finished and, with `runspec`, made by the command and from the in.gff of that spec"""
    pairs = shard_pairs(plan, e_root, outgroups, key)
    missing = [lab for _, lab in pairs if not _is_done(run_dir_of(key, lab, e_root, outgroups))]
    if missing:
        raise StopRun(f"{key}: {len(missing)} of {len(pairs)} run(s) are not finished: {missing[:6]}")
    if runspec is not None:
        foreign = [lab for run_label, lab in pairs if not _is_done(run_dir_of(key, lab, e_root, outgroups), runspec.key(key, run_label, lab))]
        if foreign:
            raise StopRun(f"{key}: {len(foreign)} run(s) were made with another command or from another in.gff than the registered one (their done.key differs): {foreign[:6]}")
    return [lab for _, lab in pairs]


def classify_group(plan: list[dict], components, e_root, outgroups, indexes, runspec: RunSpec | None = None) -> dict:
    """the call of every such record of `components`, from the finished runs of `plan` in the three outgroups (nothing else is read)"""
    known = {r["id"] for run in plan for r in run["records"]}
    lifted, unmapped = {}, {}
    for key in outgroups:
        labels = finished_shard_labels(plan, e_root, outgroups, key, runspec)
        try:
            lifted[key], unmapped[key] = X.read_lifted_runs(os.path.join(str(e_root), "lift", outgroups[key]["species"]), known, labels)
        except StopRun:
            raise
        except RuntimeError as e:                    # a finished run without its unmapped.txt, or with a table that cannot be read: damaged evidence is a stop (exit 3), not a defect
            raise StopRun(f"{key}: {e}; do not guess, report") from e
    return X.classify_components(components, lifted, unmapped, indexes)


def load_indexes(outgroups: dict = OUTGROUPS, strict: bool = False) -> dict:
    """the in_place records of each outgroup's own self-lift (read from disk)"""
    return {k: X.OutgroupIndex(X.read_selflift(t["selflift"]), k, strict_single_copy=strict) for k, t in outgroups.items()}


def _merge_anomalies(*parts: dict) -> dict:
    out: dict = {}
    for part in parts:
        for o, d in part.items():
            for k, v in d.items():
                out.setdefault(o, {})[k] = out.get(o, {}).get(k, 0) + v
    return out


def evaluate_controls(w: dict, e_root, outgroups, indexes, which=("pos", "neg"), runspec: RunSpec | None = None) -> dict:
    """Gate 7 from the finished control runs: 'pos' = the positive control (one class holds every qualifying LRPAP1 record, on the expected triple; judged on the
    multiplicity flags of the qualifying records only, X.control_view: `flagged_members` and `left_out` report the others), 'neg' = the planted negative control (at
    most one of the fake pairs is called; `records_with_triple`, `pairs_with_both_triples` and `dropped_pairs` report what it saw, without a threshold). `valid` is set
    when both were evaluated. Verdict-set runs are never read. `anomalies` counts the records of the control runs with inconsistent rows (X.anomaly_counts)."""
    out, found = {}, []
    if "pos" in which:
        calls = classify_group(w["pos"], w["control"], e_root, outgroups, indexes, runspec)
        found.append(X.anomaly_counts(calls))
        res = X.call_expansions(w["control"], calls, prefix="EC")                  # the registered rule: this is the control E-table (etable_control.tsv)
        view, flagged = X.control_view(w["ctrl"], w["control"], calls)             # the control is judged on the flags of its qualifying records only (G1)
        res_c = res if view is calls else X.call_expansions(w["control"], view, prefix="EC")
        out["positive"] = X.check_positive_control(w["ctrl"], res_c, view)
        out["positive"]["triple"] = list(next((e["triple"] for e in res_c["expansions"] if e["id"] == out["positive"]["expansion"]), ()) or [])
        out["positive"]["flagged_members"] = flagged
        out["positive"]["left_out"] = [f["record"] for f in flagged if f["left_out"]]
        out["_control_expansions"] = res["expansions"]
    if "neg" in which:
        calls = classify_group(w["ctl"], w["neg"], e_root, outgroups, indexes, runspec)
        found.append(X.anomaly_counts(calls))
        res = X.call_expansions(w["neg"], calls, prefix="EN")
        out["negative"] = X.check_negative(res, len(w["neg"]))
        out["negative"]["called_pairs"] = [dict(component=e["component"], members=[m["id"] for m in e["members"]], triple=list(e["triple"]))
                                           for e in res["expansions"]]
        out["negative"].update(X.negative_counts(w["neg"], calls, res))
    out["anomalies"] = _merge_anomalies(*found)
    if "positive" in out and "negative" in out:
        out["valid"] = bool(out["positive"]["pass"] and out["negative"]["pass"])
    return out


def control_digest(w: dict, e_root, outgroups, group: str) -> str:
    """sha256 over the outputs (done.key, loci.tsv, unmapped.txt) of every control run of `group` in the three outgroups, halves of split runs included"""
    plan = w["pos"] if group == "pos" else w["ctl"]
    h = hashlib.sha256(group.encode())
    for key in sorted(outgroups):
        for lab in shard_labels(plan, e_root, outgroups, key):
            d = run_dir_of(key, lab, e_root, outgroups)
            h.update(f"|{key}|{lab}|".encode())
            for fn in ("done.key", "loci.tsv", "unmapped.txt"):
                f = d / fn
                h.update(f.read_bytes() if f.exists() else b"-")
                h.update(b"\0")
    return h.hexdigest()


def anomaly_text(an: dict) -> str:
    return "; ".join(f"{o}: " + ", ".join(f"{k} {v}" for k, v in sorted(d.items())) for o, d in sorted(an.items()))


def flagged_text(p: dict) -> str:
    """the records of the control components that carry the multiplicity flag, with the outgroups, and which of them were left out of the control (G1)"""
    members = p.get("flagged_members") or []
    if not members:
        return ""
    text = "flagged multiplicity (record: outgroups): " + ", ".join(f"{m['record']}: {','.join(m['outgroups'])}" for m in members)
    if p.get("left_out"):
        text += "; left out of the control (not qualifying records, never a reason to fail it): " + ", ".join(p["left_out"])
    return text


def failure_text(res: dict) -> str:
    bits = []
    if "positive" in res and not res["positive"]["pass"]:
        p = res["positive"]
        if p.get("triple_ok") is False:         # one class holds every qualifying record, on the wrong loci (nothing is missing)
            msg = ("positive control (LRPAP1): the class that holds every qualifying record has the triple "
                   f"({', '.join(p.get('triple') or [])}), not the expected ({', '.join(p['expected_triple'])})")
        else:
            msg = f"positive control (LRPAP1): no single class holds every qualifying record; missing {[(m['record'], m['reason']) for m in p['missing']]}"
        flagged = flagged_text(p)
        bits.append(msg + (f"; {flagged}" if flagged else ""))
    if "negative" in res and not res["negative"]["pass"]:
        n = res["negative"]
        msg = f"negative control: {n['called']} of {n['of']} planted pairs came out as an expansion (bar {n['bar']})"
        if "records_with_triple" in n:
            msg += (f"; planted records with a triple {n['records_with_triple']}, pairs with both triples {n['pairs_with_both_triples']}, "
                    f"pairs dropped for multiplicity {n['dropped_pairs']}")
        bits.append(msg)
    if res.get("anomalies"):
        bits.append("anomalies in the Liftoff output (a record with no row that is not listed unmapped = no_row; copy rows without a primary row = no_primary_row; "
                    f"none is expected): {anomaly_text(res['anomalies'])}")
    return "; ".join(bits)


def _plausible_report(saved) -> bool:
    """a saved gate7.json worth keeping beside the new verdict: an object whose control_runs, positive and negative (when present) are objects too"""
    return isinstance(saved, dict) and all(isinstance(saved.get(k, {}), dict) for k in ("control_runs", "positive", "negative"))


def make_gate(w: dict, e_root, outgroups, load_indexes_fn, runspec: RunSpec | None = None):
    """gate(group), the function run_budgeted calls before the first unit of a group: 'pos' needs nothing, 'neg' needs the positive control to pass, 'verdict'
    needs both controls to pass. The controls are EVALUATED from the finished control runs every time the gate is asked (about 2 s): a saved verdict is never
    trusted, because a file can be written by hand. The verdict is saved in E_ROOT/gate7.json (with the digest of the control outputs it was read from) as a
    report for the reader, e.g. the triple of the positive control. Raises StopRun on a failure or on an unfinished control run."""
    status_path = Path(e_root) / "gate7.json"
    need_of = {"pos": (), "neg": ("pos",), "verdict": ("pos", "neg")}
    name_of = {"pos": "positive", "neg": "negative"}

    def gate(group: str):
        need = need_of[group]
        if not need:
            return
        try:
            saved = json.loads(status_path.read_text()) if status_path.exists() else {}
        except ValueError:
            saved = {}                                          # an unreadable report is replaced; it was never a decision
        if not _plausible_report(saved):
            saved = {}                                          # nor is one that is not even shaped like a report (a hand-edited file)
        digests = {g: control_digest(w, e_root, outgroups, g) for g in need}
        res = evaluate_controls(w, e_root, outgroups, load_indexes_fn(), which=need, runspec=runspec)
        saved.setdefault("control_runs", {})
        for g in need:
            saved["control_runs"][g] = digests[g]
            saved[name_of[g]] = res[name_of[g]]
        saved["valid"] = res.get("valid")                       # null unless both controls were evaluated in this call
        saved["anomalies"] = res["anomalies"]
        _atomic_write(status_path, json.dumps(saved, indent=1, sort_keys=True))
        if any(not res[name_of[g]]["pass"] for g in need) or res["anomalies"]:
            raise StopRun(f"Gate 7 FAILED before the {group} group: {failure_text(res)}; no verdict-set record is lifted (see {status_path})")
    return gate


# ----------------------------------------------------------------------------------------------- the E-table from the finished runs
def refuse_foreign_table(out_dir, strict: bool, hub_mode: str):
    """a table directory that already holds a table of another reading (strict single-copy or the other hub rule) is never written into"""
    sp = Path(out_dir) / "summary.json"
    if not sp.exists():
        return
    try:
        s = json.loads(sp.read_text())
    except ValueError as e:
        raise StopRun(f"{sp} is not valid JSON ({e}): choose another table directory") from e
    if s.get("strict_single_copy") != strict or s.get("hub_mode") != hub_mode:
        raise StopRun(f"{out_dir} holds a table made with strict_single_copy={s.get('strict_single_copy')} and hub_mode={s.get('hub_mode')}; this build has "
                      f"strict_single_copy={strict} and hub_mode={hub_mode}: choose another table directory")


def build_all(e_root, out_dir, *, lift_root: str | None = None, outgroups: dict = OUTGROUPS, hub_mode: str = "nodes", spec: dict | None = None, n_neg: int = N_PLANTED,
              strict: bool = False, runspec: RunSpec | None = None) -> dict:
    """Gate 7 from the finished control runs, then the E-table from the finished verdict-set runs (the control components are not in it; LRPAP1 is reported
    in gate7.json). The tables are written even when a control failed or a run has inconsistent rows; STATUS says INVALID then and `valid` is False. A control or
    verdict run that is not finished (or, with `runspec`, was made by another command or from another in.gff) is an error. A directory that holds a table of
    another reading is refused."""
    refuse_foreign_table(out_dir, strict, hub_mode)
    w = plan_from_real_inputs(hub_mode, lift_root, spec, n_neg)
    check_plan_frozen(e_root, plan_digest(w["plan"] + w["pos"] + w["ctl"], hub_mode))
    indexes = load_indexes(outgroups, strict)
    try:
        g7 = evaluate_controls(w, e_root, outgroups, indexes, runspec=runspec)
        calls = classify_group(w["plan"], w["verdict"], e_root, outgroups, indexes, runspec)
    except StopRun:
        failed = failed_controls(w, e_root, outgroups, indexes, runspec)
        if failed is None:                                     # no finished control failed: an unfinished or damaged run is a stop, as before
            raise
        return write_failed_controls(w, e_root, outgroups, out_dir, failed)
    res = X.call_expansions(w["verdict"], calls)
    anomalies = _merge_anomalies(g7["anomalies"], X.anomaly_counts(calls))
    valid = bool(g7["valid"]) and not anomalies
    gate7 = dict(valid=g7["valid"], positive=g7["positive"], negative=g7["negative"], control_components=[c.id for c in w["control"]],
                 qualifying=w["ctrl"]["qualifying"], not_required=w["ctrl"]["not_required"], planted_components=len(w["neg"]), anomalies=anomalies,
                 control_runs={g: control_digest(w, e_root, outgroups, g) for g in ("pos", "neg")})
    extra = dict(gate7_valid=g7["valid"], anomalies=anomalies, hub_mode=hub_mode, strict_single_copy=strict, verdict_components=len(w["verdict"]),
                 verdict_records=sum(len(c.such()) for c in w["verdict"]), control_components=[c.id for c in w["control"]],
                 **{f"negative_{k}": g7["negative"][k] for k in ("records_with_triple", "pairs_with_both_triples", "dropped_pairs")})
    paths = X.write_outputs(res, w["pool"], w["verdict"], calls, out_dir, rows=w["rows"], extra=extra, hub_mode=hub_mode, control=w["control"])
    X._w(os.path.join(out_dir, "etable_control.tsv"), X.ETABLE_HEADER, X.etable_rows(g7["_control_expansions"]))
    status = "VALID" if valid else "INVALID: " + failure_text(dict(g7, anomalies=anomalies))
    write_status_files(out_dir, gate7, status)
    return dict(gate7=gate7, status=status, paths=paths, expansions=len(res["expansions"]), valid=valid, dropped_classes=len(res["dropped_classes"]),
                records_in_dropped_classes=sum(len(d["members"]) for d in res["dropped_classes"]))


def write_status_files(out_dir, gate7: dict, status: str):
    """gate7.json, STATUS and SHA256SUMS (sha256sum format, every other file of the directory in name order: `sha256sum -c SHA256SUMS`)"""
    with open(os.path.join(out_dir, "gate7.json"), "w") as fh:
        json.dump(gate7, fh, indent=1, sort_keys=True)
    with open(os.path.join(out_dir, "STATUS"), "w") as fh:
        fh.write(status + "\n")
    with open(os.path.join(out_dir, "SHA256SUMS"), "w") as fh:
        for name in sorted(os.listdir(out_dir)):
            if name != "SHA256SUMS":
                fh.write(f"{X._sha(os.path.join(out_dir, name))}  {name}\n")


def failed_controls(w: dict, e_root, outgroups, indexes, runspec: RunSpec | None = None) -> dict | None:
    """The state a gate stop leaves is: the control groups lifted, the later groups never started. When a run the build needs is not finished, the finished control
    groups are evaluated one by one: the result ({positive, negative (when its runs are finished), anomalies}) if one of them FAILED or has anomalies, else None
    (then an unfinished run is a plain stop). An unfinished positive group is no evidence of anything: None."""
    try:
        pos = evaluate_controls(w, e_root, outgroups, indexes, which=("pos",), runspec=runspec)
    except StopRun:
        return None
    out, found = dict(positive=pos["positive"]), [pos["anomalies"]]
    try:
        neg = evaluate_controls(w, e_root, outgroups, indexes, which=("neg",), runspec=runspec)
        out["negative"] = neg["negative"]
        found.append(neg["anomalies"])
    except StopRun:
        pass                                                   # the planted runs were never lifted (or are damaged): nothing to say about the negative control
    out["anomalies"] = _merge_anomalies(*found)
    failed = not out["positive"]["pass"] or ("negative" in out and not out["negative"]["pass"]) or bool(out["anomalies"])
    return out if failed else None


def write_failed_controls(w: dict, e_root, outgroups, out_dir, part: dict) -> dict:
    """the table directory of a build whose control failed before the verdict runs existed: STATUS 'INVALID: <the failed control>', gate7.json and SHA256SUMS only"""
    os.makedirs(out_dir, exist_ok=True)
    groups = ("pos", "neg") if "negative" in part else ("pos",)
    gate7 = dict(valid=False, positive=part["positive"], control_components=[c.id for c in w["control"]], qualifying=w["ctrl"]["qualifying"],
                 not_required=w["ctrl"]["not_required"], planted_components=len(w["neg"]), anomalies=part["anomalies"],
                 control_runs={g: control_digest(w, e_root, outgroups, g) for g in groups})
    if "negative" in part:
        gate7["negative"] = part["negative"]
    status = "INVALID: " + failure_text(part)
    write_status_files(out_dir, gate7, status)
    return dict(gate7=gate7, status=status, paths={n: os.path.join(out_dir, n) for n in ("gate7.json", "STATUS", "SHA256SUMS")}, expansions=0, valid=False,
                dropped_classes=0, records_in_dropped_classes=0)


# ----------------------------------------------------------------------------------------------- the command line
def describe_plans(w: dict, budget_s: float) -> dict:
    """the plans and the schedule as numbers (nothing is run)"""
    units = schedule_units(w["plan"], w["pos"], w["ctl"])
    groups = group_map(w["plan"], w["pos"], w["ctl"])
    by_label = {r["label"]: r for r in w["plan"] + w["pos"] + w["ctl"]}
    pred = {u: predict_wall(len(by_label[u[1]]["records"])) for u in units}
    calls = plan_calls(units, pred, budget_s, groups)
    per_group = {}
    for g in GROUPS:
        us = [u for u in units if groups[u[1]] == g]
        per_group[g] = dict(units=len(us), predicted_work_s=round(sum(pred[u] for u in us)), calls=sum(1 for c in calls if groups[c[0][1]] == g))
    by_rtype = collections.Counter(r["rtype"] for c in w["called"] for r in c.such())
    wall = predicted_wall_by_outgroup(units, by_label)
    return dict(plan_digest=plan_digest(w["plan"] + w["pos"] + w["ctl"], w["hub_mode"]),
                called_set=dict(components=len(w["called"]), records=sum(by_rtype.values()), by_rtype=dict(by_rtype), plan=X.plan_summary(w["called_plan"])),
                control_components={c.id: len(c.such()) for c in w["control"]}, positive_plan=X.plan_summary(w["pos"]),
                planted_plan=X.plan_summary(w["ctl"]), planted_components=len(w["neg"]),
                verdict=dict(components=len(w["verdict"]), records=sum(len(c.such()) for c in w["verdict"]), plan=X.plan_summary(w["plan"])),
                units=len(units), groups=per_group, calls=len(calls), predicted_work_s=round(sum(pred.values())),
                predicted_call_s=[round(sum(pred[u] for u in c)) for c in calls], planting=X.planting_counts(w["rows"], w["pool"]),
                wall_predicted_s={o: round(v, 1) for o, v in sorted(wall.items())}, wall_factor=WALL_FACTOR,
                wall_limit_s={o: round(WALL_FACTOR * v, 1) for o, v in sorted(wall.items())},
                lrpap1=dict(qualifying=w["ctrl"]["qualifying"], not_required=w["ctrl"]["not_required"], other_non_qualifying=w["ctrl"]["other_non_qualifying"],
                            fragments=w["ctrl"]["fragments"]))


def real_root(path) -> bool:
    """is `path` the real E_ROOT, the directory of the registered E phase? A directory BELOW it (E_ROOT/rehearsal2) is a rehearsal directory."""
    return os.path.realpath(path) == os.path.realpath(E_ROOT)


def rehearsal_flags(a) -> list[str]:
    """the options that change what a run is made of (the Liftoff and minimap2 binaries, the code check, the input tables, the annotation, the reference, the size of
    the planted control, the sensitivity reading of the gate, the hub reading): legal only for a rehearsal, i.e. when --out is not the real E_ROOT"""
    bad = []
    if a.liftoff_bin != LIFTOFF_BIN:
        bad.append("--liftoff-bin")
    if a.minimap2 is not None:
        bad.append("--minimap2")
    if a.no_code_check:
        bad.append("--no-code-check")
    if a.lift_root != X.LIFT_ROOT:
        bad.append("--lift-root")
    if a.gorilla_gff != GORILLA_GFF:
        bad.append("--gorilla-gff")
    if a.gorilla_fasta != GORILLA_FASTA:
        bad.append("--gorilla-fasta")
    if a.n_planted != N_PLANTED:
        bad.append("--n-planted")
    if a.cmd == "run" and a.strict_single_copy:
        bad.append("--strict-single-copy")
    if a.rehearsal:
        bad.append("--rehearsal")
    if a.hub_mode != "nodes":
        bad.append("--hub-mode")                 # the registered rule is 'nodes'; the other reading would freeze another plan.digest on the real root
    return bad


def build_command(a) -> int:
    """`build`: Gate 7 from the finished control runs, then the E-table. Refuses (exit 2) a strict build without a table directory of its own, a table directory that
    is not empty, and the output of a rehearsal; unless --rehearsal, it verifies the run identity (code, tools, command, pinned inputs) and that every run was made by
    the registered command from its in.gff (exit 3 otherwise). Exit 4: the table was written but is INVALID."""
    deciding = os.path.join(a.out, "table")
    table_dir = a.table_dir or deciding
    if a.strict_single_copy and (not a.table_dir or os.path.realpath(table_dir) == os.path.realpath(deciding)):
        print("refusing: --strict-single-copy writes a sensitivity table and needs its own --table-dir (not <out>/table, the deciding one)", file=sys.stderr)
        return 2
    if os.path.isdir(table_dir) and os.listdir(table_dir):
        print(f"refusing: {table_dir} exists and is not empty; a build never overwrites a table (choose another --table-dir)", file=sys.stderr)
        return 2
    marker = Path(a.out) / REHEARSAL_MARKER
    if marker.exists() and not a.rehearsal:
        print(f"refusing: {a.out} holds the output of a rehearsal ({marker}); build it with --rehearsal into a table directory of its own", file=sys.stderr)
        return 2
    runspec = None
    if not a.rehearsal:
        mm2 = str(Path(a.out) / "lift" / "mm2_cache.sh")
        verify_identity(a.out, code=code_identity(), tools=tool_identity(a.liftoff_bin, MINIMAP2_ENV), command=command_fingerprint(a.liftoff_bin, a.gorilla_fasta, mm2),
                        pins=pin_inputs(input_paths(a.gorilla_fasta, a.gorilla_gff, a.lift_root, OUTGROUPS), deep=True))
        runspec = RunSpec(a.out, OUTGROUPS, a.gorilla_fasta, a.liftoff_bin, mm2)
    res = build_all(a.out, table_dir, lift_root=a.lift_root, outgroups=OUTGROUPS, hub_mode=a.hub_mode, n_neg=a.n_planted, strict=a.strict_single_copy, runspec=runspec)
    if res["valid"]:
        print(json.dumps(dict(status=res["status"], expansions=res["expansions"], dropped_classes=res["dropped_classes"],
                              records_in_dropped_classes=res["records_in_dropped_classes"], gate7=res["gate7"]["valid"]), indent=1))
        return 0
    print(json.dumps(dict(status=res["status"], gate7=res["gate7"]["valid"]), indent=1))          # no outcome statistic is printed for a table that is INVALID
    return 4


def append_run_log(path, line: str):
    """one line to E_ROOT/run.log. The file is a convenience copy of what goes to stderr, so a failure to write it is reported and not fatal (the call may have run
    Liftoff for minutes)"""
    try:
        with open(path, "a") as fh:
            fh.write(line + "\n")
    except OSError as e:
        print(f"(could not append to {path}: {e.strerror})", file=sys.stderr)


def run_command(a, w: dict, assign: dict) -> int:
    """`run`: one bounded call. Everything that can fail before a Liftoff process exists (the driver lock, the stray guard, the Liftoff binary, the reference and
    the outgroup files, the code check, the comparison with what the first call froze, the run GFFs, the links, the wrapper) comes first and writes nothing that is
    frozen; the plan and the identity are WRITTEN at the moment the first Liftoff process is about to start. The Liftoff tree inherits the lock."""
    digest = plan_digest(w["plan"] + w["pos"] + w["ctl"], a.hub_mode)
    flags = rehearsal_flags(a)
    log_path = Path(a.out) / "run.log"

    def log(msg: str):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
        print(line, file=sys.stderr)
        append_run_log(log_path, line)
    with DriverLock(a.out) as lock:
        guard_no_stray(a.out)
        check_liftoff_binary(a.liftoff_bin)
        check_supervisors()
        check_reference(a.gorilla_fasta)
        for key in sorted(OUTGROUPS):
            check_outgroup_files(OUTGROUPS[key])
        mm2_path = a.minimap2 or str(Path(a.out) / "lift" / "mm2_cache.sh")
        code = code_identity()
        tools = tool_identity(a.liftoff_bin, a.minimap2 or MINIMAP2_ENV)
        if not a.no_code_check:
            absent = [k for k, v in tools.items() if v == "absent"]
            if absent:
                raise StopRun(f"the tool(s) {', '.join(absent)} do not exist: the registered Liftoff environment is not in place")
        first = load_identity(a.out) is None
        pins = pin_inputs(input_paths(a.gorilla_fasta, a.gorilla_gff, a.lift_root, OUTGROUPS), deep=first)
        commit = a.amendment_commit                                    # a rehearsal without the code check keeps the string as typed
        if not a.no_code_check:
            needles = {**code, **{f"tool {k}": v for k, v in tools.items()}}
            if first:                                                  # the first call freezes the plan and the inputs: the committed Amendment must list THESE
                needles["plan digest"] = digest
                needles.update({f"input {n}": p["sha256"] for n, p in pins.items() if p["sha256"]})
            commit = resolve_commit(a.amendment_commit)                # 40 hex: the identity stores and compares this value, whatever the spelling of a later call
            check_committed_code(commit, files=CODE_FILES, amendment_file=AMENDMENT_FILE, needles=needles)
        extra = {"tools": tools, "command": command_fingerprint(a.liftoff_bin, a.gorilla_fasta, mm2_path)}
        freeze_identity(a.out, code, commit, write=False, extra=extra, pins=pins)
        freeze_plan(a.out, digest, write=False)
        ensure_run_gffs(a.gorilla_gff, assign, os.path.join(a.out, "gff"))
        lift_dir = make_lift_dir(a.out)
        mm2 = a.minimap2 or str(L.mm2_wrapper(lift_dir))
        units = schedule_units(w["plan"], w["pos"], w["ctl"])
        plan_by_label = {r["label"]: r for r in w["plan"] + w["pos"] + w["ctl"]}
        wall_pred = predicted_wall_by_outgroup(units, plan_by_label)
        runspec = RunSpec(a.out, OUTGROUPS, a.gorilla_fasta, a.liftoff_bin, mm2)
        gate = make_gate(w, a.out, OUTGROUPS, lambda: load_indexes(OUTGROUPS, a.strict_single_copy), runspec=runspec)

        def freeze_now():
            freeze_identity(a.out, code, commit, extra=extra, pins=pins)
            freeze_plan(a.out, digest)
            if flags:
                _atomic_write(Path(a.out) / REHEARSAL_MARKER, json.dumps({"flags": flags}, indent=1) + "\n")
        res = run_budgeted(units, plan_by_label, e_root=a.out, outgroups=OUTGROUPS, reference=a.gorilla_fasta, liftoff_bin=a.liftoff_bin, minimap2=mm2,
                           budget_s=a.budget_s, group_of=group_map(w["plan"], w["pos"], w["ctl"]), gate=gate, log=log, t0=CALL_T0, on_first_launch=freeze_now,
                           pass_fds=(lock.fd,), wall_pred_s=wall_pred)
    walls = res.get("wall_s", {})
    log("Liftoff wall so far / limit: " + ", ".join(f"{o} {walls.get(o, 0)} / {round(WALL_FACTOR * v)} s" for o, v in sorted(wall_pred.items())))
    print(json.dumps(res))
    return res["exit_code"]


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["plan", "gffs", "run", "build"])
    ap.add_argument("--hub-mode", default="nodes", choices=["nodes", "records"])
    ap.add_argument("--out", default=E_ROOT)
    ap.add_argument("--table-dir", default=None, help="build: where the E-table is written (default: <out>/table); never an existing non-empty directory")
    ap.add_argument("--amendment-commit", default=None)
    ap.add_argument("--budget-s", type=float, default=CALL_BUDGET_S)
    ap.add_argument("--strict-single-copy", action="store_true",
                    help="build only: an outgroup record that is the landing site of a copy row is not single-copy either (a sensitivity table; needs its own --table-dir)")
    ap.add_argument("--minimap2", default=None, help="rehearsals only: the minimap2 Liftoff calls (-m); fixed by Amendment 2 (default: the per-call cache wrapper around the Liftoff environment's)")
    ap.add_argument("--lift-root", default=X.LIFT_ROOT, help="the self-lift tables (gorilla/liftoff_loci.tsv, records.tsv); rehearsals only when changed")
    ap.add_argument("--gorilla-gff", default=GORILLA_GFF)
    ap.add_argument("--gorilla-fasta", default=GORILLA_FASTA)
    ap.add_argument("--liftoff-bin", default=LIFTOFF_BIN)
    ap.add_argument("--n-planted", type=int, default=N_PLANTED, help=argparse.SUPPRESS)          # registered: 100; only the tests change it
    ap.add_argument("--no-code-check", action="store_true", help="rehearsals and tests only: skip the check that the code is committed after the Amendment")
    ap.add_argument("--rehearsal", action="store_true", help="build only: the runs in --out were made by a stub (they carry REHEARSAL.json); never for the real E_ROOT")
    a = ap.parse_args(argv)
    if a.cmd == "run" and (not a.amendment_commit or not re.fullmatch(r"[0-9a-f]{7,40}", a.amendment_commit)):
        print("refusing to run Liftoff: --amendment-commit HASH of the committed Amendment 2 is required (the heavy 'go' of the user)", file=sys.stderr)
        return 2
    flags = rehearsal_flags(a)
    if a.cmd != "plan" and flags and real_root(a.out):
        print(f"refusing: {', '.join(flags)} change what a run is made of and are for rehearsals only; --out is the real root {E_ROOT} "
              f"(use another --out for a rehearsal, e.g. {E_ROOT}/rehearsal3)", file=sys.stderr)
        return 2
    if a.cmd == "build":
        return build_command(a)
    w = plan_from_real_inputs(a.hub_mode, a.lift_root, n_neg=a.n_planted)
    if a.cmd == "plan":
        print(json.dumps(describe_plans(w, a.budget_s), indent=1, sort_keys=True))
        return 0
    assign = {rc["id"]: r["label"] for r in w["plan"] + w["pos"] + w["ctl"] for rc in r["records"]}
    if a.cmd == "gffs":
        t0 = time.time()
        wrote = ensure_run_gffs(a.gorilla_gff, assign, os.path.join(a.out, "gff"))
        print(f"{'wrote' if wrote else 'kept'} {len(set(assign.values()))} run GFFs for {len(assign)} records in {time.time() - t0:.0f} s")
        return 0
    return run_command(a, w, assign)


def main(argv=None) -> int:
    """exit codes: 0 done, 75 call again, 3 a stop rule (STOP: ...), 2 refused, 4 build: the table was written but is INVALID, 70 a defect of this code (traceback),
    128 + n a signal n (the Liftoff session is killed first)"""
    try:
        return _main(argv)
    except StopRun as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 3
    except Exception:
        traceback.print_exc()
        print("INTERNAL ERROR: a defect of the driver, not a stop rule (exit 70): report it", file=sys.stderr)
        return 70


if __name__ == "__main__":
    sys.exit(main())
