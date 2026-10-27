#!/usr/bin/env python3
"""test_supplement_e2.py -- supplementary tests for the Phase E code (expansion_table.py, liftoff_cross.py), written by the second mutation review ('mut', e2).

Every test PASSES on the frozen code and was written to kill a mutant of a line that the fix round CHANGED and that the 540 tests of the fix round let survive
(see the review report: survivors and which tests kill them).  Hermetic and deterministic: made-up tables, temp directories (removed at the end of every test), stub
scripts that only `sleep`/`true`; no Liftoff, no minimap2, no real input, nothing under the real run directory is read or written (every class derives from
T.HermeticCase, which replaces E_ROOT, the self-lift root, the Liftoff binary, the gorilla FASTA and GFF, the minimap2 of the Liftoff environment and the outgroups'
files by paths that cannot exist, before any setUp runs).  No wall-clock bound anywhere: waits are polls with generous caps, never an upper bound on speed.
Only the public functions of the two modules are used.
Fix round 2 (S1: the clean-up reads /proc and signals the process group before it lists anything): the one test that pinned the old early return of kill_session was rewritten
(listed in fix_round/CHANGES.md); every other test passes as the reviewer wrote it.

Run: python3 -B -m unittest bench/hierarchy/test_supplement_e2.py     (from the repository root; it imports the fixture builders of test_expansion_table.py)
"""
import collections
import contextlib
import csv
import hashlib
import io
import json
import os
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import time
import types
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import expansion_table as X  # noqa: E402
import liftoff_cross as LC  # noqa: E402
import test_expansion_table as T  # noqa: E402  (fixture builders only)
import test_supplement_e as S  # noqa: E402  (DriverBase: the fixture of the driver tests; none of its tests is re-run here)

rec, row, lrow = T.rec, T.row, T.lrow
REAL = T.REAL


def rd(path):
    with open(path) as fh:
        return fh.read()


def write_exec(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)
    return path


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return False


class TmpCase(T.HermeticCase):
    """a temp directory per test, removed at the end of the test (the fix round's suites leave 278 of them behind per run)"""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)


# ===================================================================================================== 1. the signal handling that keeps the Liftoff tree from outliving the driver
class SignalsHeldDirectTests(TmpCase):
    """signals_held: the clean-up of a session (kill_session) is not interrupted by a second SIGTERM / SIGHUP / SIGINT. The 540 tests call it only through
    timed_watch's `finally`, and none sends a signal during the clean-up, so a mutant that holds another signal set survived all 97 tests that run it."""
    SIGS = (signal.SIGTERM, signal.SIGHUP, signal.SIGINT)

    def record(self, sig):
        seen = []
        old = signal.signal(sig, lambda s, f: seen.append(s))
        self.addCleanup(signal.signal, sig, old)
        return seen

    def test_each_of_sigterm_sighup_and_sigint_is_held_until_the_block_ends(self):
        for sig in self.SIGS:
            seen = self.record(sig)
            with LC.signals_held():
                os.kill(os.getpid(), sig)
                self.assertEqual(seen, [], signal.Signals(sig).name + " must wait for the end of the clean-up")
            self.assertEqual(seen, [sig], signal.Signals(sig).name + " is delivered when the block ends")

    def test_a_signal_that_is_not_one_of_the_three_is_not_held(self):
        seen = self.record(signal.SIGUSR2)
        with LC.signals_held():
            os.kill(os.getpid(), signal.SIGUSR2)
            self.assertEqual(seen, [signal.SIGUSR2])

    def test_the_mask_before_the_block_is_restored_and_not_cleared(self):
        before = signal.pthread_sigmask(signal.SIG_BLOCK, [signal.SIGUSR1])        # a signal the caller already blocked
        self.addCleanup(signal.pthread_sigmask, signal.SIG_SETMASK, before)
        with LC.signals_held():
            inside = signal.pthread_sigmask(signal.SIG_BLOCK, [])
            self.assertTrue({signal.SIGUSR1, *self.SIGS} <= set(inside))
        after = signal.pthread_sigmask(signal.SIG_BLOCK, [])
        self.assertIn(signal.SIGUSR1, set(after), "the caller's own mask is restored, not emptied")
        self.assertFalse(set(self.SIGS) & set(after), "the three signals are unblocked again")

    def test_the_mask_is_restored_when_the_block_raises(self):
        with self.assertRaises(RuntimeError):
            with LC.signals_held():
                raise RuntimeError("boom")
        self.assertFalse(set(self.SIGS) & set(signal.pthread_sigmask(signal.SIG_BLOCK, [])))

    def test_outside_the_main_thread_nothing_is_blocked(self):
        out = {}

        def work():
            with LC.signals_held():
                out["mask"] = set(signal.pthread_sigmask(signal.SIG_BLOCK, []))
        th = threading.Thread(target=work)
        th.start()
        th.join()
        self.assertFalse(set(self.SIGS) & out["mask"], "a worker thread is left alone (its mask is not the main thread's business)")


class TerminateOnSignalsDirectTests(TmpCase):
    SIGS = (signal.SIGTERM, signal.SIGHUP, signal.SIGINT)

    def test_each_of_the_three_signals_becomes_system_exit_128_plus_n_and_the_old_handler_comes_back(self):
        for sig in self.SIGS:
            def tripwire(signum, frame):                                                  # what would run if terminate_on_signals did not install its own handler: a failure, not the death of the test run
                raise AssertionError("terminate_on_signals did not install a handler for " + signal.Signals(signum).name)
            before = signal.signal(sig, tripwire)
            self.addCleanup(signal.signal, sig, before)
            with LC.terminate_on_signals():
                with self.assertRaises(SystemExit) as cm:
                    os.kill(os.getpid(), sig)                                             # delivered to this (single) thread before the next statement
                self.assertEqual(cm.exception.code, 128 + sig, signal.Signals(sig).name)
            self.assertIs(signal.getsignal(sig), tripwire, signal.Signals(sig).name + ": the old handler is restored")
            signal.signal(sig, before)

    def test_a_signal_that_is_not_listed_keeps_its_handler(self):
        old = signal.getsignal(signal.SIGUSR2)
        with LC.terminate_on_signals():
            self.assertIs(signal.getsignal(signal.SIGUSR2), old)

    def test_outside_the_main_thread_no_handler_is_installed_and_nothing_fails(self):
        out = {}

        def work():
            with LC.terminate_on_signals():
                out["ran"] = True
        th = threading.Thread(target=work)
        th.start()
        th.join()
        self.assertTrue(out.get("ran"))


# ===================================================================================================== 2. Gate 7 sees the anomalies of BOTH control groups
class ControlAnomalyTests(T.WorldCase):
    """evaluate_controls counts the records with inconsistent rows (no row and not listed unmapped = a lost alignment) of the positive AND of the planted runs; a mutant that
    dropped the count of the positive group survived all 73 tests that run evaluate_controls (the builder's tests make the anomaly only in verdict runs)."""

    def setUp(self):
        super().setUp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.finish_controls()

    def lose(self, label, outgroup_species, without):
        """rewrite the finished run `label` of one outgroup so that it has no row for the records in `without` and does not list them unmapped"""
        key, sp, loc = next(t for t in T.OG_SPECIES if t[1] == outgroup_species)
        run = next(r for r in self.w["pos"] + self.w["ctl"] if r["label"] == label)
        frag = {"gene-F1", "gene-F2", "gene-F3"}
        if label.startswith("pos"):
            placed = {r["id"]: (f"{loc}{2 if r['id'] in frag else 1}", 1000, 3000) for r in run["records"] if r["id"] not in without}
        else:
            j = [r["label"] for r in self.w["ctl"]].index(label)
            placed = {r["id"]: (f"{loc}{1 if j == 0 else 2}", 1000, 3000) for r in run["records"] if r["id"] not in without}
        T.make_run_dir(os.path.join(self.e_root, "lift", sp), label, placed, unmapped=[])

    def test_a_lost_record_of_the_positive_runs_is_counted(self):
        self.lose("pos2", "human", {"gene-N"})                                   # gene-N is a not-required record: the control itself still passes
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes(), which=("pos",))
        self.assertTrue(res["positive"]["pass"], res["positive"])
        self.assertEqual(res["anomalies"], {"H": {"no_row": 1}})

    def test_a_lost_record_of_the_planted_runs_is_counted(self):
        run = self.w["ctl"][0]
        self.lose("ctl1", "chimpanzee", {run["records"][0]["id"]})
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes(), which=("neg",))
        self.assertEqual(res["anomalies"], {"C": {"no_row": 1}})

    def test_both_groups_are_merged_by_outgroup_and_kind(self):
        self.lose("pos2", "human", {"gene-N"})
        self.lose("pos2", "orangutan", {"gene-N"})
        self.lose("ctl1", "human", {self.w["ctl"][0]["records"][0]["id"]})
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertEqual(res["anomalies"], {"H": {"no_row": 2}, "O": {"no_row": 1}})

    def test_the_gate_stops_on_an_anomaly_of_the_positive_group_although_the_control_passes(self):
        self.lose("pos2", "human", {"gene-N"})
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")
        self.assertIn("anomal", str(cm.exception))
        self.assertIn("H: no_row 1", str(cm.exception))

    def test_a_clean_world_has_no_anomaly_in_either_group(self):
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertEqual(res["anomalies"], {})


# ===================================================================================================== 3. the preflight of the executables
class ExecutablePreflightTests(TmpCase):
    """check_liftoff_binary / check_supervisors: the file must EXIST AND be executable (a mutant `and -> or` let a non-executable file through, so the call would freeze the
    identity and then fail at the first launch)"""

    def test_check_liftoff_binary_accepts_an_executable_file_only(self):
        good = write_exec(os.path.join(self.tmp, "liftoff"), "#!/bin/sh\n")
        LC.check_liftoff_binary(good)
        plain = os.path.join(self.tmp, "plain")
        with open(plain, "w") as fh:
            fh.write("x")
        os.chmod(plain, 0o644)
        for bad in (plain, os.path.join(self.tmp, "absent"), self.tmp):          # a file that is not executable, a missing one, a directory (executable but not a file)
            with self.assertRaises(LC.StopRun, msg=bad):
                LC.check_liftoff_binary(bad)

    def test_check_supervisors_needs_both_programs_to_exist_and_be_executable(self):
        t, tm = os.path.join(self.tmp, "timeout"), os.path.join(self.tmp, "time")
        write_exec(t, "#!/bin/sh\n")
        write_exec(tm, "#!/bin/sh\n")
        with mock.patch.object(LC, "TIMEOUT_BIN", t), mock.patch.object(LC, "TIME_BIN", tm):
            LC.check_supervisors()
        os.chmod(tm, 0o644)
        with mock.patch.object(LC, "TIMEOUT_BIN", t), mock.patch.object(LC, "TIME_BIN", tm):
            with self.assertRaises(LC.StopRun):
                LC.check_supervisors()
        os.chmod(tm, 0o755)
        os.chmod(t, 0o644)
        with mock.patch.object(LC, "TIMEOUT_BIN", t), mock.patch.object(LC, "TIME_BIN", tm):
            with self.assertRaises(LC.StopRun):
                LC.check_supervisors()
        with mock.patch.object(LC, "TIMEOUT_BIN", self.tmp), mock.patch.object(LC, "TIME_BIN", tm):
            with self.assertRaises(LC.StopRun):                                  # a directory is not a program
                LC.check_supervisors()


# ===================================================================================================== 4. the stray guard looks at PATHS, not prefixes
class StrayRootBoundaryTests(TmpCase):
    """find_stray_processes(e_root): a liftoff / minimap2 / cache wrapper of this user counts iff its working directory is the root or below it, or its command line names the
    root as a path. The six builder tests never put a process next to the root (a sibling with the same prefix), never use a symlinked root, and never test the driver's
    own process, so the mutants that turned the path rules into prefix rules or dropped one of the roots survived."""

    def setUp(self):
        super().setUp()
        self.e_root = os.path.join(self.tmp, "e")
        self.sibling = os.path.join(self.tmp, "e2")
        os.makedirs(os.path.join(self.e_root, "run"))
        os.makedirs(os.path.join(self.sibling, "run"))
        self.procs = []
        self.addCleanup(self.reap)

    def reap(self):
        for p in self.procs:
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
            p.wait()

    def spawn(self, name, cwd, *args):
        """a bash script called `name` that sleeps, started in `cwd` with `args`; returns the Popen once its command line shows the script"""
        script = os.path.join(self.tmp, "bin", name)
        os.makedirs(os.path.dirname(script), exist_ok=True)
        write_exec(script, "#!/bin/bash\nsleep 60 &\nwait\n")
        p = subprocess.Popen([script, *args], cwd=cwd, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.procs.append(p)
        t0 = time.time()
        while time.time() - t0 < 20:
            try:
                if name in rd(f"/proc/{p.pid}/cmdline"):
                    break
            except OSError:
                pass
            time.sleep(0.02)
        return p

    def pids(self, root=None):
        return [pid for pid, _ in LC.find_stray_processes(root or self.e_root)]

    def test_a_process_in_a_sibling_directory_with_the_same_prefix_is_not_ours(self):
        self.spawn("liftoff", os.path.join(self.sibling, "run"))
        self.assertEqual(self.pids(), [])

    def test_a_process_whose_cwd_is_the_root_itself_is_ours(self):
        p = self.spawn("liftoff", self.e_root)
        self.assertIn(p.pid, self.pids())

    def test_a_command_line_that_names_a_longer_path_with_the_root_as_a_prefix_is_not_ours(self):
        for other in (self.e_root + "2", self.e_root + ".bak", self.e_root + "-old", self.e_root + "_x"):
            self.spawn("minimap2", self.sibling, os.path.join(other, "genome.fa.mmi"))
        self.assertEqual(self.pids(), [])

    def test_a_command_line_that_names_the_root_as_a_path_is_ours(self):
        for named in (self.e_root, self.e_root + "/lift/human/genome.fa.mmi"):          # the root itself, a path below it (the process runs elsewhere)
            p = self.spawn("minimap2", self.sibling, named)
            self.assertIn(p.pid, self.pids(), named)

    def test_the_root_may_be_given_as_a_symbolic_link(self):
        link = os.path.join(self.tmp, "link_to_e")
        os.symlink(self.e_root, link)
        p = self.spawn("liftoff", os.path.join(self.e_root, "run"))                  # the kernel reports the resolved working directory
        self.assertIn(p.pid, self.pids(link))

    def test_the_first_word_of_the_command_line_counts_too(self):
        exe = os.path.join(self.tmp, "bin", "liftoff")
        os.makedirs(os.path.dirname(exe), exist_ok=True)
        shutil.copy("/bin/sleep", exe)
        p = subprocess.Popen([exe, "60"], cwd=os.path.join(self.e_root, "run"), start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.procs.append(p)
        t0 = time.time()
        while time.time() - t0 < 20 and "liftoff" not in rd(f"/proc/{p.pid}/cmdline"):
            time.sleep(0.02)
        self.assertIn(p.pid, self.pids())

    def test_the_process_that_asks_is_never_its_own_stray(self):
        # a child python whose own command line is `python3 .../liftoff` (a script called liftoff) with its cwd in the root asks for the strays of the root
        script = os.path.join(self.tmp, "bin", "liftoff")
        os.makedirs(os.path.dirname(script), exist_ok=True)
        with open(script, "w") as fh:
            fh.write("import os, sys\nsys.path.insert(0, %r)\nimport liftoff_cross as LC\nprint(os.getpid(), [p for p, _ in LC.find_stray_processes(%r)])\n" % (HERE, self.e_root))
        # no `-B` on the command line: the stray rule looks at the first TWO words of a command line, which must be the interpreter and the script called liftoff
        r = subprocess.run([sys.executable, script], cwd=os.path.join(self.e_root, "run"), capture_output=True, text=True, timeout=60, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        self.assertEqual(r.returncode, 0, r.stderr)
        me, found = r.stdout.split(" ", 1)
        self.assertEqual(found.strip(), "[]", f"process {me} found itself")


# ===================================================================================================== 5. the preflight of the reference and of the outgroup files: its two boundaries
class PreflightBoundaryTests(TmpCase):
    """_usable (exists AND is not empty) and the `.fai is not older than the FASTA` rule of check_reference / check_outgroup_files: boundary mutants (`> 0 -> >= 0`,
    `< -> <=`) survived the 41-48 tests that reach them, because every builder test uses a missing file or a clearly older index"""

    def genome(self, name="g.fa"):
        fasta = os.path.join(self.tmp, name)
        T.make_fake_genome(fasta)                                    # FASTA 100 s in the past, .fai and .mmi newer, k15 w10 header with one sequence per .fai line
        return fasta

    def outgroup(self, fasta):
        return dict(species="x", genome=fasta, mmi=fasta + ".mmi", selflift="/nowhere")

    def test_the_index_may_have_exactly_the_mtime_of_the_fasta_and_may_not_be_older(self):
        fasta = self.genome()
        st = os.stat(fasta)
        os.utime(fasta + ".fai", ns=(st.st_atime_ns, st.st_mtime_ns))           # the same modification time: pyfaidx keeps the index (it rebuilds only a strictly older one)
        LC.check_reference(fasta)
        LC.check_outgroup_files(self.outgroup(fasta))
        os.utime(fasta + ".fai", ns=(st.st_atime_ns, st.st_mtime_ns - 5 * 10 ** 9))
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_reference(fasta)
        self.assertIn("older", str(cm.exception))
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_outgroup_files(self.outgroup(fasta))
        self.assertIn("older", str(cm.exception))

    def test_an_empty_file_is_not_usable_whichever_of_the_three_it_is(self):
        for victim in ("", ".fai", ".mmi"):
            fasta = self.genome("v%s.fa" % victim.replace(".", "_"))
            st = os.stat(fasta)
            open(fasta + victim, "w").close()                                      # truncated to zero bytes
            os.utime(fasta, ns=(st.st_atime_ns, st.st_mtime_ns))                   # keep the order of the modification times: FASTA oldest
            os.utime(fasta + ".fai", ns=(st.st_atime_ns + 10 ** 9, st.st_mtime_ns + 10 ** 9))
            with self.assertRaises(LC.StopRun, msg="outgroup files, empty %r" % (fasta + victim)) as cm:
                LC.check_outgroup_files(self.outgroup(fasta))
            self.assertIn("missing or empty", str(cm.exception))
            if victim != ".mmi":                                                   # the reference has no index of minimap2
                with self.assertRaises(LC.StopRun, msg="reference, empty %r" % (fasta + victim)) as cm:
                    LC.check_reference(fasta)
                self.assertIn("missing or empty", str(cm.exception))

    def test_a_good_set_of_files_is_accepted(self):
        fasta = self.genome()
        LC.check_reference(fasta)
        LC.check_outgroup_files(self.outgroup(fasta))


# ===================================================================================================== 6. a launch starts from a clean run directory
class StaleFilesOfAnEarlierAttemptTests(TmpCase):
    """run_shard removes done.key, toolong, rss_stop, incomplete, loci.tsv, time.txt, lifted.gff3(.gz), unmapped.txt and in.gff_db before a launch, so that a killed or failed run
    can never be read as the run before it. The builder's test covers rss_stop, time.txt, loci.tsv and done.key; the other six names can be dropped from the list unnoticed."""
    STALE = ("done.key", "toolong", "rss_stop", "incomplete", "loci.tsv", "time.txt", "lifted.gff3", "lifted.gff3.gz", "unmapped.txt", "in.gff_db")

    def setUp(self):
        super().setUp()
        self.lifted = os.path.join(self.tmp, "stub_lifted.gff3")
        with open(self.lifted, "w") as fh:
            fh.write(T.lifted_gff([("gene-A", "h1", 1001, 2000, "+", "0.98", "0.97", 0, "False", "False", [(1001, 1400), (1701, 2000)], "gene")]))
        self.run_dir = os.path.join(self.tmp, "run")
        os.makedirs(self.run_dir)
        with open(os.path.join(self.run_dir, "in.gff"), "w") as fh:
            fh.write(T.REF_GFF)

    def make_stale(self):
        for n in self.STALE:
            with open(os.path.join(self.run_dir, n), "w") as fh:
                fh.write("wall_s\t123.0\n" if n == "time.txt" else "stale\n")

    def go(self, **kw):
        return LC.run_shard(self.run_dir, ["gene-A", "gene-B"], target="/t.fa", reference="/r.fa", minimap2="/m", env={"STUB_LIFTED": self.lifted}, **kw)

    def test_a_run_that_is_killed_at_the_hard_limit_leaves_nothing_of_the_attempt_before_it(self):
        self.make_stale()
        res = self.go(liftoff_bin=T.stub_liftoff(self.tmp, sleep=30), hard_s=1)
        self.assertEqual(res["status"], "toolong")
        self.assertIsNone(res["wall_s"], "a killed run does not report the wall of the run before it")
        left = sorted(n for n in self.STALE if os.path.exists(os.path.join(self.run_dir, n)) and rd(os.path.join(self.run_dir, n)) in ("stale\n", "wall_s\t123.0\n"))
        self.assertEqual(left, [])
        self.assertIn("killed after", rd(os.path.join(self.run_dir, "toolong")))              # the marker is this attempt's own

    def test_a_stale_toolong_or_incomplete_marker_does_not_survive_a_run_that_finishes(self):
        self.make_stale()
        res = self.go(liftoff_bin=T.stub_liftoff(self.tmp), hard_s=60)
        self.assertEqual(res["status"], "done")
        for n in ("toolong", "rss_stop", "incomplete", "lifted.gff3", "in.gff_db"):
            self.assertFalse(os.path.exists(os.path.join(self.run_dir, n)), n)
        self.assertNotEqual(rd(os.path.join(self.run_dir, "unmapped.txt")), "stale\n")
        self.assertNotEqual(rd(os.path.join(self.run_dir, "time.txt")), "wall_s\t123.0\n")


# ===================================================================================================== 7. the order of the checks before a launch
class LaunchSequenceTests(S.DriverBase):
    """run_budgeted: the stray guard runs before EVERY launch and BEFORE the freeze hook; the freeze hook runs ONCE per call, at the first launch (so a call that cannot start leaves
    nothing frozen, and a later launch does not freeze again). The builder tests only put a stray in front of the first launch."""

    def setUp(self):
        super().setUp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.units = [("H", "r01"), ("H", "r02")]

    def test_guard_then_freeze_once_then_guard_again_before_the_second_launch(self):
        events = []
        with mock.patch.object(LC, "guard_no_stray", lambda e: events.append("guard")):
            res = LC.run_budgeted(self.units, self.plan, budget_s=500, on_first_launch=lambda: events.append("freeze"), **self.kw)
        self.assertEqual(res["exit_code"], 0, res)
        self.assertEqual(events, ["guard", "freeze", "guard"])

    def test_a_stray_that_shows_up_after_the_first_run_stops_the_second_launch(self):
        calls = []

        def guard(e):
            calls.append(e)
            if len(calls) == 2:
                raise LC.StopRun("a liftoff of this root appeared")
        with mock.patch.object(LC, "guard_no_stray", guard):
            with self.assertRaises(LC.StopRun):
                LC.run_budgeted(self.units, self.plan, budget_s=500, **self.kw)
        self.assertTrue(self.done("r01"))
        self.assertFalse(self.done("r02"), "the second run must not start while a stray is alive")
        self.assertEqual(calls, [self.e_root, self.e_root])

    def test_without_the_guard_the_hook_is_still_called_once(self):
        froze = []
        res = LC.run_budgeted(self.units, self.plan, budget_s=500, stray_guard=False, on_first_launch=lambda: froze.append(1), **self.kw)
        self.assertEqual((res["exit_code"], froze), (0, [1]))

    def test_a_call_with_nothing_to_run_never_calls_the_hook_or_the_guard(self):
        LC.run_budgeted(self.units, self.plan, budget_s=500, **self.kw)                      # both runs done
        events = []
        with mock.patch.object(LC, "guard_no_stray", lambda e: events.append("guard")):
            res = LC.run_budgeted(self.units, self.plan, budget_s=500, on_first_launch=lambda: events.append("freeze"), **self.kw)
        self.assertEqual((res["exit_code"], events), (0, []))


# ===================================================================================================== 8. the kill sequence, on a simulated clock
class FakeTime:
    """LC.time replaced by a clock that moves only when the code sleeps: the kill sequence can be timed exactly, without waiting and without a wall-clock bound"""

    def __init__(self):
        self.now = 0.0

    def monotonic(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds

    strftime = staticmethod(time.strftime)
    time = staticmethod(time.time)


class KillSessionClockTests(TmpCase):
    """kill_session: SIGTERM at once, SIGKILL exactly when the grace period is over, and no waiting once the session is empty. The builder removed every wall-clock bound from the
    tests (rightly), so a mutant that always waits for the whole grace period (10 s per run, the `or -> and` of the wait loop) is caught only by a 90 s timeout."""

    def setUp(self):
        super().setUp()
        self.started = []
        self.addCleanup(self.reap)

    def reap(self):
        for p in self.started:
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
            p.wait()

    def group(self, body):
        p = subprocess.Popen(["/bin/bash", "-c", body], cwd=self.tmp, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.started.append(p)
        for _ in range(400):                                             # until the leader and its two sleeping children are up
            if len(LC.session_pids(p.pid)) >= 3:
                break
            time.sleep(0.05)
        self.assertGreaterEqual(len(LC.session_pids(p.pid)), 3)
        return p

    def kill(self, sid, **kw):
        ft, events, real_killpg = FakeTime(), [], os.killpg

        def killpg(pg, sig):
            events.append((ft.now, sig))
            return real_killpg(pg, sig)
        with mock.patch.object(LC, "time", ft), mock.patch.object(LC.os, "killpg", killpg):
            left = LC.kill_session(sid, **kw)
        return left, ft, events

    def test_a_session_that_dies_on_sigterm_is_not_waited_for_until_the_end_of_the_grace_period(self):
        p = self.group("sleep 600 & sleep 600 & wait")
        left, ft, events = self.kill(p.pid, term_grace_s=1000.0, kill_wait_s=1000.0, poll_s=100.0)
        self.assertEqual(left, [])
        self.assertEqual([sig for _, sig in events], [signal.SIGTERM])
        self.assertLessEqual(ft.now, 100.0, "at most one poll was needed; the grace period of 1000 s must not be waited out")

    def test_sigkill_comes_exactly_when_the_grace_period_is_over(self):
        p = self.group('trap "" TERM; sleep 600 & sleep 600 & wait')
        left, ft, events = self.kill(p.pid, term_grace_s=30.0, kill_wait_s=20.0, poll_s=10.0)
        self.assertEqual(left, [])
        self.assertEqual(events, [(0.0, signal.SIGTERM), (30.0, signal.SIGKILL)])
        self.assertLessEqual(ft.now, 30.0 + 20.0)

    def test_nothing_but_the_group_is_signalled_in_a_session_that_is_already_empty(self):
        # fix round 2 (S1): the process group gets its SIGTERM BEFORE anything is listed (so that a listing that fails cannot keep the group from being told). For a session that
        # is already empty that is ONE harmless SIGTERM to a group that no longer exists (ESRCH); no process is signalled, no SIGKILL follows, and nothing is waited for.
        p = subprocess.Popen(["/bin/true"], start_new_session=True)
        p.wait()
        sent, real_kill = [], os.kill
        with mock.patch.object(LC.os, "kill", lambda pid, sig: (sent.append((pid, sig)), real_kill(pid, sig))[1]):
            left, ft, events = self.kill(p.pid)
        self.assertEqual((left, events, ft.now), ([], [(0.0, signal.SIGTERM)], 0.0))
        self.assertEqual(sent, [], "no single process is signalled in an empty session")


# ===================================================================================================== 9. the run identity: every recorded component is compared
class IdentityComponentsTests(TmpCase):
    """freeze_identity / verify_identity record and compare the Amendment commit, the code hashes, the tools, the command line and the pinned inputs. The builder tests change the code,
    the tools and an input; the Amendment commit and the command line (and, for `build`, the tools and the command) have no test that changes only them."""
    CODE = {"bench/hierarchy/liftoff_cross.py": "1" * 64}
    TOOLS = {"liftoff": "a" * 64, "minimap2": "b" * 64, "gnu_time": "c" * 64}
    COMMAND = {"liftoff_bin": "/l", "reference": "/r", "minimap2": "/m", "cmd": ["x", "y"], "params": {"a": 1}}
    PINS = {"H_selflift": {"path": "/p/h", "size": 10, "mtime_ns": 5, "sha256": "d" * 64}}

    def frozen(self):
        e = os.path.join(self.tmp, "e")
        os.makedirs(e)
        LC.freeze_identity(e, self.CODE, "abc1234", extra={"tools": self.TOOLS, "command": self.COMMAND}, pins=self.PINS)
        return e

    def test_the_same_call_is_accepted_and_each_changed_component_alone_is_a_stop_that_names_it(self):
        e = self.frozen()
        same = dict(extra={"tools": self.TOOLS, "command": self.COMMAND}, pins=self.PINS)
        LC.freeze_identity(e, self.CODE, "abc1234", write=False, **same)
        changed = {
            "Amendment commit": ("abc1235", self.CODE, same),
            "code": ("abc1234", {"bench/hierarchy/liftoff_cross.py": "2" * 64}, same),
            "tools": ("abc1234", self.CODE, dict(same, extra={"tools": dict(self.TOOLS, liftoff="e" * 64), "command": self.COMMAND})),
            "command line": ("abc1234", self.CODE, dict(same, extra={"tools": self.TOOLS, "command": dict(self.COMMAND, minimap2="/other")})),
        }
        for what, (commit, code, kw) in changed.items():
            with self.assertRaises(LC.StopRun, msg=what) as cm:
                LC.freeze_identity(e, code, commit, write=False, **kw)
            self.assertIn(what, str(cm.exception))
            self.assertEqual(sum(w in str(cm.exception) for w in ("Amendment commit", "code", "tools", "command line")) >= 1, True)
        with self.assertRaises(LC.StopRun) as cm:
            LC.freeze_identity(e, self.CODE, "abc1234", write=False, extra=same["extra"], pins={"H_selflift": dict(self.PINS["H_selflift"], size=11)})
        self.assertIn("H_selflift", str(cm.exception))

    def test_build_compares_the_code_the_tools_the_command_and_the_inputs(self):
        e = self.frozen()
        good = dict(code=self.CODE, tools=self.TOOLS, command=self.COMMAND, pins=self.PINS)
        LC.verify_identity(e, **good)
        for what, kw in (("code", dict(good, code={"bench/hierarchy/liftoff_cross.py": "3" * 64})), ("tools", dict(good, tools=dict(self.TOOLS, gnu_time="f" * 64))),
                         ("command line", dict(good, command=dict(self.COMMAND, reference="/elsewhere"))),
                         ("H_selflift", dict(good, pins={"H_selflift": dict(self.PINS["H_selflift"], mtime_ns=6)}))):
            with self.assertRaises(LC.StopRun, msg=what) as cm:
                LC.verify_identity(e, **kw)
            self.assertIn(what, str(cm.exception))

    def test_the_first_call_writes_nothing_when_asked_only_to_compare_and_records_everything_otherwise(self):
        e = os.path.join(self.tmp, "e2")
        os.makedirs(e)
        LC.freeze_identity(e, self.CODE, "abc1234", write=False, extra={"tools": self.TOOLS}, pins=self.PINS)
        self.assertFalse(os.path.exists(os.path.join(e, "run_identity.json")))
        LC.freeze_identity(e, self.CODE, "abc1234", extra={"tools": self.TOOLS, "command": self.COMMAND}, pins=self.PINS)
        rec_ = json.loads(rd(os.path.join(e, "run_identity.json")))
        self.assertEqual(sorted(rec_), ["amendment_commit", "code_sha256", "command", "inputs", "tools"])
        self.assertEqual((rec_["amendment_commit"], rec_["code_sha256"], rec_["tools"], rec_["command"], rec_["inputs"]), ("abc1234", self.CODE, self.TOOLS, self.COMMAND, self.PINS))


# ===================================================================================================== 10. who counts as a member of the session
class SessionMembersTests(TmpCase):
    """session_pids lists the LIVE members of a session (a zombie is not one), and kill_session signals every member whatever its process group: Liftoff's tree can contain
    processes in other process groups of the same session (a shell with job control puts every background job in a group of its own)"""

    def setUp(self):
        super().setUp()
        self.started = []
        self.addCleanup(self.reap)

    def reap(self):
        for p in self.started:
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
            p.wait()

    def start(self, body):
        p = subprocess.Popen(["/bin/bash", "-c", body], cwd=self.tmp, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.started.append(p)
        return p

    def ps_all(self, sid):
        out = subprocess.run(["ps", "-o", "pid=,stat=", "--sid", str(sid)], capture_output=True, text=True).stdout
        return [(int(l.split()[0]), l.split()[1]) for l in out.splitlines() if len(l.split()) >= 2]

    def test_a_zombie_is_not_a_live_member_of_the_session(self):
        # a python parent that starts a child, never waits for it and then becomes `sleep` (exec): the child stays a zombie of that process for as long as it lives (bash would reap it
        # from its SIGCHLD handler, which makes a bash stand-in a race)
        script = write_exec(os.path.join(self.tmp, "zombie_parent.py"), "import os, subprocess\nsubprocess.Popen(['/bin/true'])\nos.execvp('sleep', ['sleep', '600'])\n")
        p = subprocess.Popen([sys.executable, "-B", script], cwd=self.tmp, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.started.append(p)
        zombies = []
        for _ in range(400):
            zombies = [pid for pid, st in self.ps_all(p.pid) if st.startswith("Z")]
            if zombies:
                break
            time.sleep(0.05)
        self.assertTrue(zombies, "the stand-in must have left a zombie for this test to mean anything")
        self.assertEqual(LC.session_pids(p.pid), [p.pid])
        left = LC.kill_session(p.pid, term_grace_s=1.0, kill_wait_s=1.0, poll_s=0.05)
        self.assertEqual(left, [], "a zombie is not a process that survived the kill")

    def test_members_in_other_process_groups_of_the_session_are_signalled_too(self):
        p = self.start("set -m; sleep 600 & sleep 600 & wait")           # job control: each background job gets a process group of its own, in the same session
        pids = []
        for _ in range(400):
            pids = LC.session_pids(p.pid)
            if len(pids) >= 3:
                break
            time.sleep(0.05)
        self.assertGreaterEqual(len(pids), 3, pids)
        groups = {os.getpgid(pid) for pid in pids}
        self.assertGreaterEqual(len(groups), 2, "the stand-in must have members in more than one process group")
        left = LC.kill_session(p.pid, term_grace_s=5.0, kill_wait_s=5.0, poll_s=0.05)
        self.assertEqual(left, [])
        for pid in pids:
            self.assertFalse(alive(pid), pid)


# ===================================================================================================== 11. which inputs are hashed: the two size limits are exclusive
class PinBoundaryTests(TmpCase):
    """pin_inputs: a file BELOW INPUT_BIG_BYTES carries its sha256 (when deep, else only below SHALLOW_SHA_BYTES); a file of exactly that size does not. The boundary mutant
    (`< -> <=`) survived the 43 tests that run pin_inputs; the limits are module constants, so the boundary is tested with small ones."""

    def test_a_file_of_exactly_the_big_limit_is_not_hashed_and_one_byte_less_is(self):
        for size, hashed in ((9, True), (10, False), (11, False)):
            path = os.path.join(self.tmp, "g%d.bin" % size)
            with open(path, "wb") as fh:
                fh.write(b"x" * size)
            with mock.patch.object(LC, "INPUT_BIG_BYTES", 10), mock.patch.object(LC, "SHALLOW_SHA_BYTES", 10 ** 6):
                got = LC.pin_inputs({"f": path}, deep=True)["f"]
            self.assertEqual((got["size"], got["sha256"] is not None), (size, hashed), size)
            if hashed:
                self.assertEqual(got["sha256"], hashlib.sha256(b"x" * size).hexdigest())

    def test_a_shallow_pin_hashes_only_below_the_shallow_limit_and_a_file_of_exactly_that_size_is_not_hashed(self):
        for size, hashed in ((4, True), (5, False), (6, False)):
            path = os.path.join(self.tmp, "h%d.bin" % size)
            with open(path, "wb") as fh:
                fh.write(b"y" * size)
            with mock.patch.object(LC, "INPUT_BIG_BYTES", 10 ** 6), mock.patch.object(LC, "SHALLOW_SHA_BYTES", 5):
                shallow = LC.pin_inputs({"f": path}, deep=False)["f"]
                deep = LC.pin_inputs({"f": path}, deep=True)["f"]
            self.assertEqual(shallow["sha256"] is not None, hashed, size)
            self.assertIsNotNone(deep["sha256"], "a deep pin hashes every file below the big limit")


# ===================================================================================================== 12. the registered hub count of summary.json: "two such records" includes exactly two
class HubCountBoundaryTests(TmpCase):
    """summary.json's `hubs_with_two_such_records` is the registered 7 (hubs among the components with AT LEAST two placed records of 200 bp or more). The builder's test has hubs with
    three and with one such record, never with exactly two, so `>= 2 -> > 2` (a registered number off by every hub of two records) survived."""

    def summary(self, hub_ids):
        rows, recs = T.mini_selflift()
        pool = X.build_pool(rows, recs)
        for cid in hub_ids:
            pool.components[cid].sites += ["copy@x%s:%d-%d+" % (cid, i, i + 9) for i in range(60)]          # more than 51 nodes: a hub
        called = X.called_components(pool, "nodes")
        calls = {r["id"]: dict(triple=T.T1, reason=None, per={k: dict(n_place=1, detail="") for k in "HCO"}) for c in called for r in c.such()}
        out = os.path.join(self.tmp, "out")
        X.write_outputs(X.call_expansions(called, calls), pool, called, calls, out, hub_mode="nodes")
        return pool, json.loads(rd(os.path.join(out, "summary.json")))

    def test_a_hub_with_exactly_two_such_records_is_counted(self):
        pool, s = self.summary(["PC_gene-D"])                                      # gene-D and gene-E: two placed records of 2 kb
        self.assertEqual(len(pool.components["PC_gene-D"].such()), 2)
        self.assertEqual((s["hubs"], s["hubs_with_two_such_records"]), (1, 1))

    def test_a_hub_with_one_such_record_is_a_hub_but_not_counted_and_one_with_three_is(self):
        pool, s = self.summary(["PC_gene-F", "PC_gene-A"])
        self.assertEqual((len(pool.components["PC_gene-F"].such()), len(pool.components["PC_gene-A"].such())), (1, 3))
        self.assertEqual((s["hubs"], s["hubs_with_two_such_records"]), (2, 1))

    def test_a_component_of_exactly_the_hub_size_is_not_a_hub_and_one_more_node_is(self):
        rows, recs = T.mini_selflift()
        pool = X.build_pool(rows, recs)
        c = pool.components["PC_gene-D"]
        c.sites += ["copy@z:%d-%d+" % (i, i + 9) for i in range(X.HUB_MAX - c.size("nodes"))]
        self.assertEqual(c.size("nodes"), X.HUB_MAX)
        self.assertFalse(c.hub("nodes"))
        c.sites.append("copy@z:999-1008+")
        self.assertTrue(c.hub("nodes"))


# ===================================================================================================== 13. the clean-up of a session is shielded from a second signal (end to end)
SHIELD_STANDIN = """#!/bin/bash
echo "$$ $(ps -o pgid= -p $$ | tr -d ' ') $(ps -o sid= -p $$ | tr -d ' ')" > "$PWD/pids.txt"
trap 'echo termed >> "$PWD/termed.txt"' TERM
sleep 600 &
sleep 600 &
while :; do sleep 0.1; done
"""


class CleanUpShieldTests(TmpCase):
    """timed_watch kills the Liftoff session in a `finally` under signals_held. A stand-in whose leader traps SIGTERM (and survives it until SIGKILL) keeps the clean-up busy for the
    whole grace period; a SECOND SIGTERM to the driver in that window must wait for the end of the clean-up. Without the shield (`with signals_held():` removed) the second signal
    raises SystemExit inside kill_session, the driver exits at once and the session is left alive: the 8 to 14 GB tree that outlives the driver (the documented crash).
    No test of the builder or of the mutation review sends a second signal."""

    def setUp(self):
        super().setUp()
        self.procs = []
        self.addCleanup(self.reap)

    def reap(self):
        for p in self.procs:
            if p.poll() is None:
                p.kill()
            p.communicate()
        for pid in os.listdir("/proc"):                                    # whatever the stand-in left in a session started under this test's temp directory
            if pid.isdigit():
                try:
                    if os.readlink(f"/proc/{pid}/cwd").startswith(self.tmp):
                        os.kill(int(pid), signal.SIGKILL)
                except (OSError, ProcessLookupError):
                    pass

    def test_a_second_sigterm_during_the_clean_up_waits_for_its_end(self):
        stand_in = write_exec(os.path.join(self.tmp, "standin.sh"), SHIELD_STANDIN)
        driver = write_exec(os.path.join(self.tmp, "driver.py"), T.DRIVER_SCRIPT)
        e_root = os.path.join(self.tmp, "e")
        run_dir = os.path.join(e_root, "run")
        p = subprocess.Popen([sys.executable, "-B", driver, HERE, e_root, stand_in, "60", "5", "3", "nolock"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.procs.append(p)
        pids_file = os.path.join(run_dir, "pids.txt")
        for _ in range(600):                                                # until the stand-in has recorded its ids and its two sleepers are up
            if os.path.exists(pids_file) and rd(pids_file).split():
                sid = int(rd(pids_file).split()[2])
                if len(LC.session_pids(sid)) >= 3:
                    break
            time.sleep(0.05)
        else:
            self.fail("the stand-in did not start")
        p.send_signal(signal.SIGTERM)                                       # the first signal: the handler raises SystemExit and the clean-up begins
        termed = os.path.join(run_dir, "termed.txt")
        for _ in range(600):                                                # the stand-in's leader has seen the SIGTERM of kill_session: the clean-up is under way (3 s of grace)
            if os.path.exists(termed):
                break
            time.sleep(0.05)
        else:
            self.fail("the clean-up did not start")
        p.send_signal(signal.SIGTERM)                                       # the second signal, inside the grace period
        out, err = p.communicate(timeout=120)
        self.assertEqual(p.returncode, 128 + signal.SIGTERM, err)
        self.assertEqual(LC.session_pids(sid), [], "the session must be empty when the driver has exited")


# ===================================================================================================== 14. the numbers the plan command prints
class PlanNumbersTests(T.WorldCase):
    """describe_plans prints the predicted Liftoff wall of every outgroup and its limit to ONE decimal (the Amendment quotes 1,669.1 s and 3,338.2 s). The builder's test uses
    figures that are whole numbers, so rounding to no decimal (`round(v, 1) -> round(v, 0)`) survived."""

    def setUp(self):
        super().setUp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_the_predicted_wall_and_the_limit_are_printed_with_one_decimal(self):
        d = LC.describe_plans(self.w, 530)
        units = LC.schedule_units(self.w["plan"], self.w["pos"], self.w["ctl"])
        by_label = {r["label"]: r for r in self.w["plan"] + self.w["pos"] + self.w["ctl"]}
        want = {o: sum(LC.predict_wall(len(by_label[lab]["records"])) for oo, lab in units if oo == o) for o in "HCO"}
        self.assertNotEqual({o: round(v, 1) for o, v in want.items()}, {o: round(v) for o, v in want.items()}, "the fixture must discriminate one decimal from none")
        self.assertEqual(d["wall_predicted_s"], {o: round(v, 1) for o, v in sorted(want.items())})
        self.assertEqual(d["wall_limit_s"], {o: round(LC.WALL_FACTOR * v, 1) for o, v in sorted(want.items())})
        self.assertEqual(d["wall_factor"], LC.WALL_FACTOR)


# ===================================================================================================== 15. a record without a placement has no multiplicity (E-table side)
class NoPlacementIsNoMultiplicityTests(T.HermeticCase):
    """has_multiplicity: `a record that has no placement has no multiplicity`. classify_record counts the copy rows of a record even when its primary row is not a placement, so only the
    early return of has_multiplicity keeps such a record out of the dropped classes, out of record_calls.tsv's multiplicity column and out of the control pool (as `outgroup_copy_row`,
    eligible 1, instead of its own reason `no_placement`, eligible 0). The mutant `return False -> return True` survived the 144 tests that call it."""

    def setUp(self):
        self.idx = T.make_indexes()
        rows, recs = T.mini_selflift()
        self.pool = X.build_pool(rows, recs)
        self.comp = self.pool.components["PC_gene-D"]                          # gene-D, gene-E: two placed records of 2 kb
        self.unmapped = {"H": set(), "C": set(), "O": set()}

    def calls(self):
        lifted_e = {k: {"gene-E": v} for k, v in T.lifted_for("gene-E", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)).items()}
        lifted_d = {k: {"gene-D": v} for k, v in T.lifted_for("gene-D", ("h2", 1000, 3000), ("c2", 1000, 3000), ("o2", 1000, 3000)).items()}
        # gene-D in H: the primary row is a weak partial mapping (no placement) and the record also has a reverse copy row of coverage 0.9
        lifted_d["H"]["gene-D"] = [T.lrow("gene-D", 0, "h2", 1000, 3000, cov=0.4, sid=0.9, partial=1), T.lrow("gene-D_1", 1, "x9", 0, 800, cov=0.9, sid=0.9)]
        return {"gene-D": X.classify_record("gene-D", lifted_d, self.unmapped, self.idx), "gene-E": X.classify_record("gene-E", lifted_e, self.unmapped, self.idx)}

    def test_the_record_has_no_placement_in_h_and_no_multiplicity_there(self):
        d = self.calls()["gene-D"]
        h = d["per"]["H"]
        self.assertEqual((h["status"], h["copy_rows"], h["multiplicity"]), ("no_placement", 1, False))
        self.assertIs(X.has_multiplicity(h), False)
        self.assertEqual((d["multiplicity"], d["reason"], d["triple"]), ([], "no_placement", None))

    def test_has_multiplicity_on_hand_made_calls(self):
        for info, want in ((dict(status="no_placement", copy_rows=3, n_place=2, multiplicity=True), False), (dict(status="no_placement", detail="no_primary_row", copy_rows=1, n_place=1), False),
                           (dict(status="ok", copy_rows=1, n_place=1), True), (dict(status="ok", copy_rows=0, n_place=2), True), (dict(status="ok", copy_rows=0, n_place=1, multiplicity=True), True),
                           (dict(status="ok", copy_rows=0, n_place=1, multiplicity=False), False), (dict(status="outgroup_copy_row", copy_rows=0, n_place=1), False), ({}, False)):
            self.assertIs(X.has_multiplicity(info), want, info)

    def test_such_a_record_neither_drops_a_class_nor_enters_the_control_pool_nor_the_multiplicity_column(self):
        calls = self.calls()
        res = X.call_expansions([self.comp], calls)
        self.assertEqual((res["expansions"], res["dropped_classes"]), ([], []))
        self.assertEqual(X.control_pool_rows([self.comp], calls, res, {}), [], "gene-D has no placement in H: it is not flagged, so the component has no control-pool reason")
        out = os.path.join(tempfile.mkdtemp(), "out")
        self.addCleanup(shutil.rmtree, os.path.dirname(out), True)
        X.write_outputs(res, self.pool, [self.comp], calls, out, hub_mode="nodes")
        with open(os.path.join(out, "record_calls.tsv")) as fh:
            table = list(csv.reader(fh, delimiter="\t"))
        by_gene = {r[1]: r for r in table[1:]}
        self.assertEqual(by_gene["gene-D"][-1], "", "the multiplicity column of a record without a placement is empty")
        self.assertEqual(json.loads(rd(os.path.join(out, "summary.json")))["records_with_reverse_copy_rows"], {"H": 0, "C": 0, "O": 0})

    def test_the_same_copy_row_on_a_record_that_is_placed_does_flag_it(self):
        lifted_d = {k: {"gene-D": v} for k, v in T.lifted_for("gene-D", ("h2", 1000, 3000), ("c2", 1000, 3000), ("o2", 1000, 3000)).items()}
        lifted_d["H"]["gene-D"].append(T.lrow("gene-D_1", 1, "x9", 0, 800, cov=0.9, sid=0.9))
        d = X.classify_record("gene-D", lifted_d, self.unmapped, self.idx)
        self.assertEqual((d["per"]["H"]["status"], d["multiplicity"]), ("ok", ["H"]))
        self.assertIs(X.has_multiplicity(d["per"]["H"]), True)


# ===================================================================================================== 16. a session that cannot be emptied is a stop
class UnkillableSessionTests(TmpCase):
    """timed_watch raises StopRun when kill_session reports a process that survived SIGKILL (`do not start another run until ps --sid is empty`). No test reaches that branch:
    the mutant `if left: raise StopRun -> if False:` survived the 99 tests that run timed_watch."""

    def test_a_process_that_survives_sigkill_is_a_stop_that_names_the_session_and_the_pid(self):
        stand_in = write_exec(os.path.join(self.tmp, "ok.sh"), "#!/bin/bash\nexit 0\n")
        d = os.path.join(self.tmp, "run")
        os.makedirs(d)
        with mock.patch.object(LC, "kill_session", lambda sid, **k: [424242]):
            with self.assertRaises(LC.StopRun) as cm:
                LC.timed_watch([stand_in], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=30, poll_s=0.1)
        self.assertIn("424242", str(cm.exception))
        self.assertIn("not empty after SIGKILL", str(cm.exception))

    def test_an_empty_session_after_the_kill_is_the_normal_case(self):
        stand_in = write_exec(os.path.join(self.tmp, "ok.sh"), "#!/bin/bash\nexit 0\n")
        d = os.path.join(self.tmp, "run")
        os.makedirs(d)
        rc, peak = LC.timed_watch([stand_in], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=30, poll_s=0.1)
        self.assertEqual(rc, 0)


# ===================================================================================================== 17. "more than one placement": exactly two rows of coverage >= 0.5 are enough
class PlacementCountBoundaryTests(T.HermeticCase):
    """classify_record flags a record `multiplicity` in an outgroup when it has a row with extra_copy_number > 0 OR MORE THAN ONE row of coverage >= 0.5 (registered reading 5). Real Liftoff
    output always has the first condition when it has the second, so no builder test gives a record two primary-tagged rows; the boundary `n_place > 1 -> n_place > 2` survived."""

    def setUp(self):
        self.idx = T.make_indexes()
        self.unmapped = {"H": set(), "C": set(), "O": set()}

    def call(self, h_rows):
        lifted = {k: {"g1": v} for k, v in T.lifted_for("g1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)).items()}
        lifted["H"]["g1"] = h_rows
        return X.classify_record("g1", lifted, self.unmapped, self.idx)

    def test_two_rows_of_coverage_one_half_flag_the_record_and_keep_the_locus_of_the_primary_one(self):
        c = self.call([T.lrow("g1", 0, "h1", 1000, 3000), T.lrow("g1", 0, "h2", 1000, 3000, cov="0.5")])
        h = c["per"]["H"]
        self.assertEqual((h["n_place"], h["copy_rows"], h["status"], h["locus"], h["multiplicity"], h["detail"]), (2, 0, "ok", "gene-H1", True, "reverse_copy_row"))
        self.assertEqual((c["multiplicity"], c["triple"]), (["H"], ("gene-H1", "gene-C1", "gene-O1")))
        self.assertTrue(X.has_multiplicity(h))

    def test_a_second_row_just_below_one_half_does_not_count(self):
        h = self.call([T.lrow("g1", 0, "h1", 1000, 3000), T.lrow("g1", 0, "h2", 1000, 3000, cov="0.49")])["per"]["H"]
        self.assertEqual((h["n_place"], h["multiplicity"], X.has_multiplicity(h)), (1, False, False))

    def test_one_row_is_no_multiplicity_and_three_rows_are(self):
        h = self.call([T.lrow("g1", 0, "h1", 1000, 3000)])["per"]["H"]
        self.assertEqual((h["n_place"], h["multiplicity"]), (1, False))
        h = self.call([T.lrow("g1", 0, "h1", 1000, 3000), T.lrow("g1", 0, "h2", 1000, 3000), T.lrow("g1", 0, "h3", 1000, 3000)])["per"]["H"]
        self.assertEqual((h["n_place"], h["multiplicity"]), (3, True))


if __name__ == "__main__":
    unittest.main()
