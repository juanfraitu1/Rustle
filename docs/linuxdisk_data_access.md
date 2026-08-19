# Linux data-disk access

The large reference assemblies and alignment data live on a separate 1.8 TB ext4 disk. They are
not stored in the Rustle repository and must never be deleted or modified as part of a benchmark
cleanup.

## Mount or reattach after a WSL restart

The disk is Windows `PHYSICALDRIVE0`, partition 2. `wsl --mount` does not survive a WSL restart.
The commands below are copied from the guarded mount instructions in `/home/juanfra/.bashrc`.
The attach/detach syntax is also documented by
[Microsoft's WSL disk-mount guide](https://learn.microsoft.com/en-us/windows/wsl/wsl2-mount-disk).

First, in an **Administrator PowerShell** window:

```powershell
wsl --mount \\.\PHYSICALDRIVE0 --partition 2
```

Then, inside WSL:

```bash
sudo mkdir -p /mnt/linuxdisk
sudo mount --bind /mnt/wsl/PHYSICALDRIVE0p2 /mnt/linuxdisk
```

Do not hard-code `/dev/sde2` or another `/dev/sd*` name: the kernel device letter can change between
sessions. Use the stable WSL-generated `PHYSICALDRIVE0p2` path.

Verify the mount before reading data or starting a large build:

```bash
mountpoint -q /mnt/linuxdisk
findmnt -T /mnt/linuxdisk
test -r /mnt/linuxdisk/home/juanfraitu/gorilla_haps/mat.fa
```

The shell startup guard sets `CARGO_TARGET_DIR` to
`/mnt/linuxdisk/home/juanfraitu/rustle_target` only while the mount is active. If the disk is absent,
the variable is deliberately unset so a large Rust build does not silently fill the Windows C:
drive.

In a managed or sandboxed session the mount may be exposed read-only even when it is normally
writable. Check `findmnt -T /mnt/linuxdisk` and `test -w /mnt/linuxdisk` before selecting an output
directory. Never attempt to remount it from an automated task. Read reference data in place and
write small durable certificates under the Rustle `bench/` tree; use `/tmp` only for disposable
intermediates.

## Ape and human assembly inventory

The following files were verified readable on 2026-08-17:

| species/haplotype | assembly file | index/status |
|---|---|---|
| human CHM13v2.0 | `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa` | reference used by the HSA O1 audit |
| gorilla KB3781 maternal | `/mnt/linuxdisk/home/juanfraitu/gorilla_haps/mat.fa` | `mat.fa.fai`; 225 sequences |
| gorilla KB3781 paternal | `/mnt/linuxdisk/home/juanfraitu/gorilla_haps/pat.fa` | `pat.fa.fai`; 24 sequences |
| chimpanzee AG18354 | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna` | `.fai`; 26 sequences |
| Bornean orangutan AG05252 | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna` | `.fai`; 26 sequences |

The corresponding optional annotation files are `GGO_genomic.gff`, `PTR_genomic.gff`, and
`PPY_genomic.gff` under `winloci_data/`. They may be used to label results after inference, but the
O1 rooting rules must not require them.

Existing gorilla minimap2 indexes include:

```text
/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.mat.asm20.rebuild.mmi
/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.pat.asm20.mmi
/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.mat.splice.mmi
/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.pat.splice.mmi
```

The `asm20` indexes are suitable for sensitive duplicated-block discovery. Strict `asm5` flank
mapping should use `mat.fa`/`pat.fa` directly or a separately generated `asm5` index: minimap2 index
seed parameters are fixed when the index is built.

## Detach safely

After all processes using the disk have stopped, unmount the bind mount in WSL:

```bash
sudo umount /mnt/linuxdisk
```

Then detach it in Administrator PowerShell:

```powershell
wsl --unmount \\.\PHYSICALDRIVE0
```
