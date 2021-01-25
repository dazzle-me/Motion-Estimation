# Motion-Estimation

Motion estimation is the algorithm of lossy video sequence comperession which works relying on the fact that two adjacent frames are quite similar, and it is pretty easy to reconstruct consecutive frame from the previous one.

Here is the result of algorithm execution, 

on the top - original video sequence, 

on the bottom - every ```i-th``` frame was derived as some reconstruction from the ```(i - 1)-th``` frame from the original video. 

![](motion_estimation.gif)

## Build
```
python setup.py build_ext -i
```
## Algorithm 
Follow ``` motions_estimation.ipynb```
