To contribute code to pyjetty, you should make a pull request from your fork of the repo.

For example, working with a single branch:

1. First, create your fork of `matplo/pyjetty` by clicking the "Fork" tab in the upper right of the `matplo/pyjetty` webpage.

2. Next, go to your local pyjetty repository (doing `git clone` from the original `matplo/pyjetty`). We want to add a "remote" to your fork. To see your current remotes, you can type `git remote -v`. 

- To add a remote to your fork, type `git remote add myFork git@github.com:<github-username>/pyjetty.git` (or using https, if you use that)

3. Suppose you now have a commit on your local repository that you want to propose to add to the main `pyjetty` repository.

- First, we should update our local code to the latest changes in `matplo/pyjetty` and put our commit on top of that: `git pull --rebase`
- Next, commit to your fork: `git push myFork master`
- Then, go to your fork and create a pull request.
