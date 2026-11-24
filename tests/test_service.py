import unittest
from wavelength_calibration import process


class CalibrationTests(unittest.TestCase):
    def base(self):
        return {"max_residual_nm":"0","reference_lines":[{"id":"l0","pixel":"0","wavelength_nm":"400"},{"id":"l1","pixel":"10","wavelength_nm":"421"},{"id":"l2","pixel":"20","wavelength_nm":"444"},{"id":"l3","pixel":"30","wavelength_nm":"469"}],"query_pixels":["5","15"]}

    def test_recovers_quadratic(self):
        result=process(self.base())
        self.assertEqual(result["coefficients"],{"c0":"400","c1":"2","c2":"1/100"})
        self.assertEqual(result["calibrated"][0]["wavelength_nm"],"1641/4")

    def test_zero_residual_is_accepted(self):
        self.assertEqual(process(self.base())["status"],"accepted")

    def test_residual_can_reject_fit(self):
        data=self.base(); data["reference_lines"][2]["wavelength_nm"]="445"
        self.assertEqual(process(data)["status"],"rejected")

    def test_query_outside_range_fails(self):
        data=self.base(); data["query_pixels"]=["31"]
        with self.assertRaises(ValueError): process(data)

    def test_lines_must_be_monotone(self):
        data=self.base(); data["reference_lines"][2]["wavelength_nm"]="420"
        with self.assertRaises(ValueError): process(data)

    def test_singular_pixels_fail(self):
        data=self.base(); data["reference_lines"][1]["pixel"]="0"
        with self.assertRaises(ValueError): process(data)


if __name__=="__main__": unittest.main()

