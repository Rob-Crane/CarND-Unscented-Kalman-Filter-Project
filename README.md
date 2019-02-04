# Unscented Kalman Filter Project Starter Code

Project submission for Udacity Term 2 Unscented Kalman Filter.

## Evaluating Consistency

Setting the environment variable `NIS_DATA_DIR` will cause normalized innovation squared (NIS) values to be saved to a comma-separated values (CSV) file in the specified directory.  The name of the file will encode the current, hard-coded values for longitudunal accelaration and yaw accelaration.  A MATLAB script, included in project root folder, will generate summary plots and table for every NIS CSV file in directory passed as argument to `plotNIS` function.
