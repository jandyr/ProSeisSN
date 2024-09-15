----
# Processamento do Sinal Sísmico (ProSeisSN)

This repository contains relevant material and notebooks to guide students through the practical activities of the Course **ProSeisSN**.

----
## Running the course codes on the cloud
- **Run** the notebooks from this **GitHub** Repo on the cloud by opening **Binder** on another tab \
  by clicking on $\rightarrow$ 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jandyr/ProSeisSN/HEAD)

- **Share** the **Binder** with the URL:\
     `https://mybinder.org/v2/gh/jandyr/ProSeisSN/HEAD`
----

## ProSeisSN structure

- **Unit**$\mathbf{nn}$: Material for each Unit of the Course \
 **├── MdlExpo**: Expositive material in *pdf* \
 **└──  MdlNbk**:  Jupyter Notebooks used in the practical units \
 **$\quad$ $\quad$ ├── Data**: Data used in the Notebooks \
 **$\quad$ $\quad$ └── Codes**: Extra codes for the Notebooks 
---
##---------------------- Repository  -------------------------------
 
## Course Environment

- **JupyterLab**: a web--based user interface, which works with {\bf Jupyter Notebooks}, containing computer code, and rich text elements.

- **Binder**: to open and execute the {\bf Jupyter Notebooks} on the cloud.

- **Google Colab**: a better alternative to run the notebooks from github.

- **GitHub**: a *Git* cloud repository with a version control system. You can clone the *Repo* to your computer.

- **GitHub Classroom**: for storing your assignments and Projects under your course’s {\bf GitHub} course organization.

- **Discord**: for discussions and questions about the lessons and exercises.

- **Anaconda**: to create a Python environment on your machine for running the Course's codes locally.

## Course content

- [**Unit01**](https://github.com/jandyr/ProSeisSN/tree/main/Unit01): Overview of Geophysics. Course environment.

- [**Unit02**](https://github.com/jandyr/ProSeisSN/tree/main/Unit02): Introduction to Obspy and time series analysis.



- **Notebook handling**: You can ue **Colab** by going to  [Google Colab](https://colab.research.google.com/), which supports URLs that link directly to a GitHub browser.
 
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

Students are directed to reference from the Literature in the **Unit Notes**. The following refers to the basics in contructing the Course environment.

* [Short Python tutorial](https://swcarpentry.github.io/python-novice-inflammation/index.html)

* [Longer Python tutorial](https://docs.python.org/3/tutorial/index.html)
* [Getting Started with Anaconda](https://docs.anaconda.com/anaconda/user-guide/getting-started/)
* [Jupyter Notebook Overview](https://jupyter-notebook.readthedocs.io/en/stable/)
* [ObsPy Tutorial](https://docs.obspy.org/tutorial/)
