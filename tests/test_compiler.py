import unittest
import tempfile
import time
import os
from pathlib import Path

from moments2amplitudes.compiler import (
    ProcessLock,
    get_aclic_target_so_paths,
    get_aclic_source_paths,
    is_aclic_updated,
    clean_aclic_artifacts,
    ensure_aclic_updated,
    TARGET_LIBRARIES,
)


class TestCompiler(unittest.TestCase):
    def test_target_libraries_list(self):
        self.assertIn("ParameterHelper", TARGET_LIBRARIES)
        self.assertIn("Equation", TARGET_LIBRARIES)
        self.assertIn("EquationSolver", TARGET_LIBRARIES)
        self.assertIn("MomentHelper", TARGET_LIBRARIES)
        self.assertIn("MassDependentFitter", TARGET_LIBRARIES)

    def test_process_lock_acquisition_and_release(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            lock_path = Path(tmp_dir) / "test.lock"

            # Primary acquisition
            with ProcessLock(lock_path, timeout=1.0):
                self.assertTrue(lock_path.is_file())

                # Secondary acquisition from same or other process should time out if non-blocking
                lock2 = ProcessLock(lock_path, timeout=0.2, wait_interval=0.05)
                with self.assertRaises(TimeoutError):
                    with lock2:
                        pass

            # Once exited, another lock can be acquired immediately
            with ProcessLock(lock_path, timeout=1.0):
                self.assertTrue(lock_path.is_file())

    def test_is_aclic_updated_missing_so(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            src_dir = root / "src"
            src_dir.mkdir(parents=True)
            # Create a source file
            (src_dir / "ParameterHelper.cpp").write_text("// test")

            updated, reason = is_aclic_updated(root)
            self.assertFalse(updated)
            self.assertIn("Missing target library", reason)

    def test_is_aclic_updated_outdated_source(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            src_dir = root / "src"
            src_dir.mkdir(parents=True)

            # Create all target .so files with an older timestamp
            past_time = time.time() - 100
            for lib in TARGET_LIBRARIES:
                p = src_dir / f"{lib}_cpp.so"
                p.write_bytes(b"dummy binary")
                os.utime(p, (past_time, past_time))

            # Create a source file with newer timestamp
            src_file = src_dir / "ParameterHelper.cpp"
            src_file.write_text("// updated source")
            current_time = time.time()
            os.utime(src_file, (current_time, current_time))

            updated, reason = is_aclic_updated(root)
            self.assertFalse(updated)
            self.assertIn("modified after", reason)

    def test_is_aclic_updated_success(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            src_dir = root / "src"
            src_dir.mkdir(parents=True)

            past_time = time.time() - 100
            now_time = time.time()

            # Create source files with older timestamp
            src_file = src_dir / "ParameterHelper.cpp"
            src_file.write_text("// source")
            os.utime(src_file, (past_time, past_time))

            # Create all target .so files with newer timestamp
            for lib in TARGET_LIBRARIES:
                p = src_dir / f"{lib}_cpp.so"
                p.write_bytes(b"dummy binary")
                os.utime(p, (now_time, now_time))

            updated, reason = is_aclic_updated(root)
            self.assertTrue(updated)
            self.assertIn("up-to-date", reason)

    def test_clean_aclic_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            src_dir = root / "src"
            src_dir.mkdir(parents=True)

            # Create dummy artifacts
            (src_dir / "ParameterHelper_cpp.so").write_bytes(b"so")
            (src_dir / "ParameterHelper_cpp.d").write_text("d")
            (src_dir / "ParameterHelper_cpp_ACLiC_dict_rdict.pcm").write_bytes(b"pcm")
            (src_dir / "ParameterHelper_cpp_ACLiC_dict.cxx").write_text("cxx")
            (src_dir / "ParameterHelper_cpp_tmp_123").write_text("tmp")
            # Create a persistent source file that should NOT be cleaned
            (src_dir / "ParameterHelper.cpp").write_text("// cpp")

            removed = clean_aclic_artifacts(root)
            self.assertEqual(len(removed), 5)
            self.assertFalse((src_dir / "ParameterHelper_cpp.so").exists())
            self.assertFalse((src_dir / "ParameterHelper_cpp.d").exists())
            self.assertFalse((src_dir / "ParameterHelper_cpp_ACLiC_dict.cxx").exists())
            self.assertTrue((src_dir / "ParameterHelper.cpp").exists())

    def test_ensure_aclic_updated_fast_path(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            src_dir = root / "src"
            src_dir.mkdir(parents=True)

            now_time = time.time()
            for lib in TARGET_LIBRARIES:
                p = src_dir / f"{lib}_cpp.so"
                p.write_bytes(b"dummy")
                os.utime(p, (now_time, now_time))

            # If libraries are up-to-date, ensure_aclic_updated returns True without invoking root
            res = ensure_aclic_updated(repo_root=root)
            self.assertTrue(res)


if __name__ == "__main__":
    unittest.main()
