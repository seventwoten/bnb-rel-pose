# bnb-rel-pose

bnb-rel-pose implements branch-and-bound search for estimating relative pose without input correspondences in Matlab.
Currently handles three cases:
- Rotation only with no translation
- Translation with known rotation
- Rotation and translation (test does not search full range of R and t yet)

## Test scripts

Run the test scripts in Matlab, with all other .m files in the same folder:
- test_R_only.m
- test_t_only.m
- test_Rt.m
