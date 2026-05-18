import unittest
from unittest.mock import patch

from msspack.doctor import run_doctor


class DoctorTests(unittest.TestCase):
    def test_run_doctor_marks_ume_optional(self) -> None:
        with patch("msspack.doctor.which", side_effect=["/usr/bin/java", None]), patch(
            "msspack.doctor._importable",
            return_value=True,
        ), patch("msspack.doctor.list_installed", return_value={}):
            checks = run_doctor()

        by_name = {check.name: check for check in checks}
        self.assertFalse(by_name["BUSCO (optional)"].ok)
        self.assertFalse(by_name["DDBJ parser"].ok)
        self.assertFalse(by_name["DDBJ transchecker"].ok)
        self.assertTrue(by_name["DDBJ ume (optional)"].ok)
