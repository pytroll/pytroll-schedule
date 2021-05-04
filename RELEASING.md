# Releasing pytroll-schedule

1. checkout main branch
2. pull from repo
3. run the unittests
4. run `loghub` and update the `CHANGELOG.md` file:

```
loghub pytroll/pytroll-schedule -u <username> -st v0.5.0 -plg bug "Bugs fixed" -plg enhancement "Features added" -plg documentation "Documentation changes"
```

Don't forget to commit!


5. Create a tag with the new version number, starting with a 'v', eg:

```
git tag v0.5.1 -m "Version 0.5.1"
```

See [semver.org](http://semver.org/) on how to write a version number.

6. push changes to github `git push --follow-tags`

7. Verify travis tests passed and deployed sdist and wheel to PyPI
