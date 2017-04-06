START_DIR="${PWD}"
START_DIR=`echo $START_DIR | sed 's/\ /\\ /g'`

# Binds to the project directory when using bash
PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR=`echo $PROJECT_DIR | sed 's/\ /\\ /g'`

echo "Appending project directory to PATH variable."
PATH=$PROJECT_DIR:$PATH

# Adds a line to the bash profile that always adds the project directory
# to the $PATH environment variable
if [ -f "$HOME/.bash_profile" ] ; then
  echo "

# Added by build-construct setup script
PATH=\"${PROJECT_DIR}:\$PATH\"" >> $HOME/.bash_profile

elif [ -f "$HOME/.profile" ] ; then
  echo "

# Added by build-construct setup script
PATH=\"${PROJECT_DIR}:\$PATH\"" >> $HOME/.profile

else
  echo "

# ~/.bash_profile: executed by the command interpreter for login shells.

# Added by build-construct setup script
PATH=\"${PROJECT_DIR}:\$PATH\"" >> $HOME/.bash_profile
fi



cd $PROJECT_DIR

# Clones the custom logger class into the project directory
if ! [ -d "${PROJECT_DIR}/log" ] ; then
  echo "Cloning logging repo."
  git clone "https://github.com/njherrmann/log.git" log
else
  echo "Logging repo found."
fi


# Adds log to the python path variable on startup
echo "Adding logging repo to python path."
echo "

# Added by build-construct setup script
export PYTHONPATH=\$PYTHONPATH:${PROJECT_DIR}/log" >> $HOME/.bashrc

# Reloads the .bashrc file (with added PYTHONPATH entry)
source ~/.bashrc



# easy_installs requests, bs4, lxml on OS X systems
# pip installs the same modules on Linux systems
echo "Installing python modules."
if [[ $OSTYPE == darwin* ]] ; then
  sudo easy_install requests bs4 lxml
elif [[ $OSTYPE == linux* ]] ; then
  pip install requests bs4 lxml
fi


# Creates the pair-guides script that simply runs the python script.
echo "Building master script."
echo $PROJECT_DIR/src/pair-guides.py | sed 's/\ /\\ /g' > "$PROJECT_DIR/pair-guides"
chmod +x "$PROJECT_DIR/pair-guides"


cd "$START_DIR"