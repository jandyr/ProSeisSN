
# Processamento do Sinal Sísmico (ProSeisSN)

This repository contains relevant material and notebooks to guide students through the practical activities of the Course **ProSeisSN**.

## Repository Organization

The repository is organized as follows.

- **Data**: input data used in the Notebooks.

- **Notebooks**: Python codes and Jupyter Notebooks used in the practical units.

- **Unit Notes**: Written material for each Unit of the Course.

## Course Environment

- **JupyterLab**: a web--based user interface, which works with {\bf Jupyter Notebooks}, containing computer code, and rich text elements.

- **Binder**: to open and execute the {\bf Jupyter Notebooks} on the cloud. 

- **GitHub**: a *Git* cloud repository with a version control system. You can clone the *Repo* to your computer.

- **GitHub Classroom**: for storing your assignments and Projects under your course’s {\bf GitHub} course organization.

- **Discord**: for discussions and questions about the lessons and exercises.

- **Anaconda**: to create a Python environment on your machine for running the Course's codes locally.

## Course content

- [**Unit01**](https://github.com/jandyr/ProSeisSN/tree/main/Unit01): Overview of Geophysics. Course environment.

- [**Unit02**](https://github.com/jandyr/ProSeisSN/tree/main/Unit02): Introduction to Obspy and time series analysis.


## Running the course codes

- **On the cloud**: Use *Binder* by clicking on 

 [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jandyr/ProSeisSN/main/)

To test, click on

     `hello.ipynb`

 to run the run {\bf Jupyter Notebook} on the cloud.
 
- **Locally**: Create the *Course Environment* on your local machine, by installing the relevant software. Follow the instructions in the *Unit Notes* to install

* Anaconda or Miniconda;
* GitHub;
* Jupyter Notebook and JupyterLab;
* Python, ObsPy etc.

Fork this repository to a directory on your local machine by

     `git clone https://github.com/jandyr/ProSeisSN`
     `cd ProSeisSN`

this way your local changes are kept *local*. Alternatively you can simply *Download ZIP* to your local directory, using the *big green button*.

Create a conda local environment for running data exercises using

     `conda env create -f environment.yml`.

Activate the environment with

     `conda activate ProSeisSN`

and deactivate it with

     `conda deactivate`

## Additional Resources

Students are directed to reference from the Literature in the **Unit Notes**. The following refer to the basic to contruct the Course environment.

* [Short Python tutorial](https://swcarpentry.github.io/python-novice-inflammation/index.html)

* [Longer Python tutorial](https://docs.python.org/3/tutorial/index.html)
* [Getting Started with Anaconda](https://docs.anaconda.com/anaconda/user-guide/getting-started/)
* [Jupyter Notebook Overview](https://jupyter-notebook.readthedocs.io/en/stable/)
* [ObsPy Tutorial](https://docs.obspy.org/tutorial/)
