#! /bin/bash

if [ $# -lt 2 ]
then
    GRAPHIZ_DOT_FILES=log/*.gv
    GRAPHIZ_PS_FILES_PATH=log
else
    GRAPHIZ_DOT_FILES=$@
    GRAPHIZ_PS_FILES_PATH=$(dirname $1)
fi

if [ ! -e $GRAPHIZ_DOT_FILES ]
then
    echo "there's no file in there..."
    exit 1
fi


for file in $GRAPHIZ_DOT_FILES
do
    GVFILE=$(basename $file .${file##*.})
#    if [ -f $GRAPHIZ_PS_FILES_PATH/$GVFILE".png" ];
#    then
#        echo "File $GVFILE.png exists"
#    else
#        echo "render file" $GVFILE".png"
#        dot -Tpng $file -o $GRAPHIZ_PS_FILES_PATH/$GVFILE".png"
#    fi
    if [ -f $GRAPHIZ_PS_FILES_PATH/$GVFILE".ps" ];
    then
      echo "File $GVFILE.ps exists"
    else
      echo "render file" $GVFILE".ps"
      dot -Tps $file -o $GRAPHIZ_PS_FILES_PATH/$GVFILE".ps"
    fi
done
