# stremr 1.1.2 (2026-02-13)

## New Features
* Added new routine `get_RD_diff()` to formally test effect modification
* Extended the output of `get_RDs()` so it can return the influence curve values to enable the comparison of two risk differences between two baseline subgroups (effect modification test)
* Generalized analytic estimates of IPW and TMLE variance to allow clustering of longitudinal observations (e.g., multiple cohort entries from the same patient in trial emulation)
* Revised the estimator of the IPW standard error so that it returns inflated values with near positivity violation when using saturated MSM. Original implementation continues to be used if the MSM is not saturated or for the unadjusted/crude estimator.
* Revised the estimator of the TMLE standard error so that it returns inflated values with near positivity violation when user requests to cap residuals so they can never be smaller than a user-specified (or a default) lower bound. This avoids under-inflated standard error estimates in situations where extreme large weights are multiplied by extremely small Q residuals.

## Bug Fixes
* Fixed a bug where stabilized inverse probability weights could be NaN instead of 0 when the numerator was 0
* Fixed various minor code errors identified when extending the scope of applications to several real-world comparative effectiveness studies

## Minor Changes
* All influence curve values are now stamped with the cluster ID

# stremr 1.1.3 (2026-02-16)

## Bug Fixes
* Ensured that the revised estimator of the IPW standard error was only used when calling survMSM with a saturated MSM

## Minor Changes
* Updated `README` and `.Rbuildignore`
* Added `Makefile` to automate rendering of both HTML and Markdown documentation
* Added `NEWS` to track version changes over time
