# DLBclass Classification Tool

This tool, written in Python, consists of a lightweight module in the form of a 
Jupyter notebook for visually aided classification. There are two major steps that need to be followed:

1) Setup - this includes cloning the repository, virtual environment setup, and package installation.

2) Running the Jupyter notebook: DLBclass_classify.ipynb

### Setup

1) Users must have `Python3` and `git` installed on their machine. If you already have
   both `git` and `Python3` installed, move to the step #2 within the setup.
    
    <br>

    <b>Installing `git`</b>
    
    open a terminal/command prompt and enter the command `git --version`. 
    If your OS says something like "command not found", git needs to be installed,
    and so follow the instructions for your OS here:  

   https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

   <br>
   
   <b>Installing `Python3`</b>
   
   Next, in the terminal enter `python --version`, and ensure that your version of Python is at least Python 3.7.x.
   
   If not, install the latest python version at:
   
   https://www.python.org/downloads/

   <br>

2) Clone the repository to your local machine by entering the following commands:

   <br>
   <b>Mac / Linux users:</b>

   `cd ~/Desktop`
   
   `git clone git@github.com:getzlab/DLBclass-tool.git`
   
   `cd DLBclass-tool`
   
   <br>
   
   <b>Windows users:</b>
   
   `chdir Desktop`
   
   `git clone git@github.com:getzlab/DLBclass-tool.git`
   
   `chdir DLBclass-tool`

   <br>
3) Setup your virtual environment. Enter the command:

   `python -m venv ./venv`

   <br>
   
4) Activate your virtual environment. Enter the command:

    <b>Mac / Linux users:</b>
    
    `source ./venv/bin/activate`
    
    <b>Windows users:</b>
    
    `.\venv\Scripts\activate.bat`

   <br>
   
5) Install required packages. This may take up to 5 minutes or so.

   `pip install -r requirements.txt`

### Running the classification notebook

   In the same command line session, enter the command:

   `jupyter notebook`

   once your browser opens the directory, click on the notebook named

   `DLBclass_classify.ipynb`

   <b>Follow the instructions within the notebook to classify your samples.</b>
   
   



