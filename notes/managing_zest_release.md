# Zest installer

## First login

1. Make sure you login using git commit -am "message" and git push to check

2. Version of number in setup.py must be the most recent.

3. Package actually gets released in the github action test, see github action for more information.

4. Delete obsolete wheels, might disturb the release.

5. Note that both requirements.txt should be updated with the right dist as well as
install_requires in setup.py. Make sure to update them both with the same dependencies.