# Synopsis

[pyjetty](https://github.com/matplo/pyjetty) (currently private) uses functionality implemented in [heppy](https://github.com/matplo/heppy) (public)

# Contributing

- follow the pull request scheme - pushing to branch master must be reviewed - if you are unsure what's that please ask - see more explanation below

Now that we are starting to share a base of common code -- and since we are soon producing publication results with it -- we should be a bit more prudent about our git workflow. It's important that anyone using the shared classes can expect that there won't be any changes to the shared code without their knowledge.

This concerns the following code locations:
alice_analysis/analysis/user/substructure
alice_analysis/process/base
alice_analysis/analysis/base
I propose that we move to a Pull Request workflow for these shared directories. To be explicit, this means that instead of committing directly to master, you should do the following:
Push your commit(s) to a branch other than master on the remote. For example, if I have my commits in a local branch my-dev, I should do: git push origin my-dev:my-dev
Make a pull request from my-dev to master through github's web interface. You should tag anyone who might be affected, and add any relevant comments. 
Allow some time for others to comment / code review before merging the PR to master.
I also noticed that you guys are sometimes doing merge commits. I realize there is some legitimate philosophical differences on when to use merge commits, but I would like to ask that we avoid doing merge commits for this shared code wherever possible -- these make it more difficult to understand code changes (for example the github diff on PRs does not handle these well). Our developments are sufficiently slow that merge commits are anyway really not necessary -- I strongly suggest we use rebase instead, where necessary. 
