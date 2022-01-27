# Debugging s.hgt

## Step 1:

Generate pvals plot for one of the body parts to find where the s.hgt is NA.

## Step 2:

Starting from fit object, run arbutus step by step to see where the NA starts.

## Step 3:

Identify how to stop NA from being NA. Maybe default to some value?

Matts Note: Add error to the fitContinuous function.

Note: Adding error doesn't seem to fix the problem.
