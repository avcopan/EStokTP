# run with: . /path/to/fake-install.sh
export THIS_DIR="`pwd`/`dirname "$BASH_SOURCE"`"
export PROJECT_ROOT=$THIS_DIR/..

export PATH=$THIS_DIR/bin:$PROJECT_ROOT/exe:$PATH
