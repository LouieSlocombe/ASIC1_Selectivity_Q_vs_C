#!/bin/bash

git add *.py
git add *.pdf
git add */*.pdf
git add */*/*.pdf
git add */*/scan.pickle
git add */*/*.html

git commit -m "auto-commit"
git push
