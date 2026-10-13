# Construction fiducials

One `.mrk.json` per subject, holding the five fiducials `Nz`, `Iz`, `LPA`, `RPA`
and `Cz` from which rule `landmarks` builds the full 10-10 system on the
subject's scalp (`workflow/scripts/make_landmarks.py`). The sample sheet's
`landmarks` column points at the file; `Cz` is optional, as the landmark builder
re-estimates the cranial vertex anyway.

Coordinates are read in the **RAS of the subject's staged NIfTI**, in mm — the
space cedalion reads `headmodel/landmarks.mrk.json` in. The `coordinateSystem`
field has to say so: `RAS` (or an `MNI` variant for a template whose RAS is an
MNI space). A file exported from 3D Slicer with its default `LPS` is converted;
anything else is rejected rather than guessed at.

For your own subject, pick the four fiducials in a viewer and export them. Note
that enabling ACPC preprocessing changes the staged image's RAS, so a file
picked on the original image no longer applies; the rule errors out when a
fiducial lies further than `freesurfer.landmarks.max_snap_mm` from the scalp,
which is what that mistake looks like.

## The two files shipped here

Only the five fiducials are kept, not the full landmark sets they come from.

| file | source |
| --- | --- |
| `colin27.mrk.json` | FieldTrip's `standard_1005.elc`, which is defined on the Colin27 head. Taken from cedalion's copy (`cedalion.data.get_fieldtrip_colin27_landmarks`). Colin27 is distributed registered to MNI305, hence `coordinateSystem: mni305`. |
| `icbm152.mrk.json` | Cutini, Moro and Bisconti (2010), *Functional near infrared optical imaging in cognitive neuroscience: an introductory review*, Table 1. doi [10.1016/j.neuroimage.2010.09.030](https://doi.org/10.1016/j.neuroimage.2010.09.030). The ICBM152 templates are in MNI152, hence `coordinateSystem: mni152`. |

The ICBM152 values were transcribed from the table as it is reproduced in
cedalion's `examples/head_models/48_headmodel_landmarks_verification.ipynb`, not
from the paper itself.

Both sets are also what cedalion's own published head models were built from, so
the fiducials here match theirs; the landmarks built from them do not, because
they are constructed on this pipeline's FreeSurfer scalp rather than on
cedalion's SPM-derived one.
