# vg documentation

## Method write-ups

Long-form descriptions of methods live here in the source tree, beside the code that implements
them, so that a method and its description are reviewed in one diff and versioned together. A
write-up covers what a user needs to understand and configure a method: the objective it
optimises, every term and parameter, and what its outputs do and do not mean. Rationale that only
matters to someone *editing* the code stays in the headers.

- [read-likelihood-genotyping.md](read-likelihood-genotyping.md) — the `vg call --read-likelihood`
  model: objective function, all terms, all parameters, the linkage layer, and the VCF fields.

The [wiki](https://github.com/vgteam/vg/wiki) remains the home for tutorials and worked examples.
It is a separate repository (mounted here as the `wiki` submodule), so anything published there
cannot be reviewed alongside a code change; a page there should link to the write-up rather than
restate it.

# Automated markdown manpage

Make a markdown document with the usage messages of (selected) `vg` subcommands. 
Calls the `vg` command, so it will match the version available in the command line.

```sh
python3 vgmanmd.py > wiki/vg-manpage.md
```

Then commit and push the changes to the wiki submodule, or copy the markdown content to the [wiki page](https://github.com/vgteam/vg/wiki/vg-manpage).

## Edit descriptions

See [`vgmanmd.desc.md`](vgmanmd.desc.md) file. 
This file must be in the current directory. 
Also, in the title `# <NAME>`, `<NAME>` should match name of the command. 
The name `intro` is used for the introduction, to be placed after the table of content.
The names don't need to match a command, and not all commands have to be described: if available, a subcommand description will be added before its usage message.

The goal of the description is to be vague enough that we don't need to change them often, but informative enough that new users can get a good feel of the commands, plus pointers to other resources (e.g. Wiki pages).

## Change list of selected subcommands

At the top of [`vgmanmd.py`](vgmanmd.py), change the *cmds* list.

