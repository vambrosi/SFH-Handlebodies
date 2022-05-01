## SFH of Handlebodies

This is a python program that computes the graded sutured Floer homology of a sutured handlebody.

This is a work in progress. So far the code computes SFH for a handlebody of genus 2 in the case where the train track giving the sutures is generic (triangle inequality is satisfied).

For examples, and more details see the Jupyter notebook [HandlebodyExamples](./HandlebodyExamples.ipynb).

The computation is done by working with a fixed maximal family of compressing disks and then computing the bordered sutured Floer homology of each of the resulting pieces (all 3-balls with parametrized disks on their boundaries). Some of the code in the modules.py file is due to Jonathan Hanselman.
