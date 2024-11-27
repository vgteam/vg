# Automated markdown manpage

Make a markdown document with the usage messages of (selected) `vg` subcommands. 
Calls the `vg` command, so it will match the version available in the command line.

```sh
python3 vgmanmd.py > man.md
```

Then copy the markdown content to a wiki page, for example, or move the markdown page to somewhere in the vg repo?

## Edit descriptions

See [`vgmanmd.desc.md`](vgmanmd.desc.md) file. 
This file must be in the current directory. 
Also, in the title `# <NAME>`, `<NAME>` should match name of the command. 
The name `intro` is used for the introduction, to be placed after the table of content.
The names don't need to match a command, and not all commands have to be described: if available, a subcommand description will be added before its usage message.

The goal of the description is to be vague enough that we don't need to change them often, but informative enough that new users can get a good feel of the commands, plus pointers to other resources (e.g. Wiki pages).

## Change list of selected subcommands

At the top of [`vgmanmd.py`](vgmanmd.py), change the *cmds* list.

