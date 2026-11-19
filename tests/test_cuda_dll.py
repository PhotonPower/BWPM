"""Windows: directories searched for the CUDA DLLs (cuFFT, cudart) before importing bwpm._core.

On Windows the RPATH of the Linux build does not exist, and Python ≥ 3.8 does not search PATH for
the dependencies of extension modules. `bwpm.fast` therefore adds the CUDA directories of the
build (written by CMake), or CUDA_PATH if the build recorded none.
"""

from __future__ import annotations

from pathlib import Path

from bwpm.fast import _cuda_dll_directories


def test_build_directories_win_over_cuda_path(tmp_path: Path) -> None:
    built = tmp_path / "v13.5" / "bin"
    (built / "x64").mkdir(parents=True)
    env = tmp_path / "v13.4"
    (env / "bin" / "x64").mkdir(parents=True)
    dirs = _cuda_dll_directories([str(built), str(tmp_path / "gone")], str(env))
    assert dirs == [str(built), str(built / "x64")]


def test_cuda_path_when_the_build_recorded_none(tmp_path: Path) -> None:
    (tmp_path / "bin").mkdir()
    dirs = _cuda_dll_directories([str(tmp_path / "gone")], str(tmp_path))
    assert dirs == [str(tmp_path / "bin")]


def test_nothing_without_build_and_environment() -> None:
    assert _cuda_dll_directories([], None) == []
