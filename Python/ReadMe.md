## Requirements
Make sure that you are using Python version 3.11.4.
Make sure that the necessary libraries are installed.
These libraries can be installed and then the script, with the proper imports, can be run in Python.
There are multiple ways to install these libraries:
	-Use Python's direct installation package installer "pip"
		EX: pip install numpy
	-Use Ananconda
		EX: conda install numpy"
	-Use virtual environments. Specifically Python's venv module which can be used to create isolate Python environments.
		EX: 
		python -m venv myenv
		source myenv/bin/activate
		pip install numpy
	-There are a variety of other methods that work as well. Feel free to use whichever bets suits your needs.
Make sure to import the proper functions from their libraries and the proper functions from the other Python files. Many of the functions written in the python files are called within other python files and must be imported first.

## Possible Errors
If the code runs but the output does not reflect the imbedded output images:
	-Make sure the variables are calculated correctly and there are no small typos(typos should not be present to start). Check the calculations for variables in other scripts to ensure they match.
	-Ensure the 'Params['NeuronPopulation']' parameter is properly set to the amount of neurons for the given code --> this applies to only the WTA and synfire scripts.
		-4 for WTA
		-3 for synfire
	-Check '# Run the actual sim'
        offset --> This should be set to 10
		Params['Input'][:, 0] --> This should be


