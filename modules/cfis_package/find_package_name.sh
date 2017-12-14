echo "Below is a list of every file in the template where package_name is mentioned:"
echo
grep -r "package_name" . --exclude="*build*" --exclude="*.pyc"  --exclude="*template*"
