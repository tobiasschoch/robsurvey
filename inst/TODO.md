# MILESTONE R 0.2
## ADD FEATURE
* vignette: robust regression
* model-assisted estimation
    - variance estiamtor
    - vignette
* tests
* svymean_reg_huber (+ _tukey, + _total_)
    - functionality
    - types:
        - projective
        - ADU
        - psi
        - Lee
        - Duchesne
    - summary + variance
* .cov_reg_compund: regression: variance: compound design-model distribution
* plot.svyreg_rob
* svymean_: examples add:
    - scale()
    - summary()
    - coef()
    - SE()
    - residuals()
* svyreg
    - predict

## FIX
* makefile
* survey weighted regression: model-based cov: BUG

## DOCUMENTATION
* svymean_huber and weighted_mean_huber (+ _tukey, + _total_)
    - RHT vs Hajek
    - rwm vs. rht => discuss in weighted_mean_huber (+ _tukey); set a pointer
      from svymean_huber
* svymean_reg + svymean_reg_huber + svymean_reg_tukey (+ _total_)
    - description
    - details
    - examples
* svyreg
    - description
    - details
    - examples
* svyreg_huber (+ _tukey)
    - details
    - examples

# MILESTONE R 0.3

## DOCUMENTATION
* weighted_line + weighted_median_lines
    - check if "cars" data example is usable


# UNSCHEDULED

## NEW FEATURES
* robust calibration
* Dalen survey method

