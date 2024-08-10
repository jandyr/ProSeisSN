## Course-coppe

# Processamento do Sinal e do Ruído Sísmicos: Aplicações à Engenharia e ao Meio--Ambiente (ProSeisSN)

## Content

The repository is organized as follows:

- **Data**: input data used in the Notebooks.

- **Notebooks**: Python codes and Jupyter Notebooks used in the practical units:

- [**Obspy**](https://github.com/https://github.com/jandyr/ProSeisSN/main/Obspy/ObspyIntro.ipynb): Introduction to Obspy and time series analysis.


## Environment

The `environment.yml` file 
To ensure reproducibility of the results, we have provided an . Ensure to have installed Anaconda or Miniconda on your computer. 



You can view and run the course materials online by clicking on the badge below:

Launch it here: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jandyr/ProSeisSN/main/)

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/jandyr/ProSeisSN/main

## Course content
Unit-01 - Blabla

## Installation and usage in your local machine
1) Install Anaconda or Miniconda on your computer.

2) Using github:
2.1. Fork this repository to your local machine by
     `git clone url`; the url is given by the *big green button*

1. If u
(forking will allow you to keep your local changes seperate from changes here);
    - Clone ()
1. If just downloading:
    - Download using the *big green button*

1) Use
```
conda env create -f environment.yml
```
to create a Python environment using conda for running data exercises.

2) To activate the environment,
```
conda activate ProSeisSN
```
and deactivate it with
```
conda deactivate
```



    
## How to use these notebooks

2. Change into the newly created directory;
3. Install the requirements (recommended to use [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#id2)):
    ```bash
    conda env create -f environment.yml  # Create an environment called gphs445
    source activate gphs445  # Activate the environment - You will need to do this everytime you use the notebooks
    ```
4. Start jupyter-lab (run `jupyter lab`) in your repository directory and navigate to the notebook you 
   want to work on;
5. Save any changes you make in the notebook;
6. When you want to back-up your changes, or when you are happy with them, commit the
   changes and push to your upstream repository 
   (check out the [github cheatsheet](https://services.github.com/on-demand/downloads/github-git-cheat-sheet.pdf) for more commands):
   ```bash
   git add <notebook you have been working on> # Replace with the filename you were working on
   git commit -m "Some memorable commit message"
   git push origin master
   ```




## Additional Resources
* [Python tutorial](https://docs.python.org/3/tutorial/index.html)
* [Getting Started with Anaconda](https://docs.anaconda.com/anaconda/user-guide/getting-started/)
* [Jupyter Notebook Overview](https://jupyter-notebook.readthedocs.io/en/stable/)
* [ObsPy Tutorial](https://docs.obspy.org/tutorial/)


## Contact

Please direct any questions on  to [wavefrontgeo@gmail.com](mailto:wavefrontgeo@gmail.com).

