### Install YouCompleteMe

It seems that in China, the installation process by Vundle is blocked. I got around it by manual git. This works on os x and arch linux.

```bash
cd ~/.vim/bundle
git clone git@github.com:Valloric/YouCompleteMe.git
cd YouCompleteMe
git submodule update --init --recursive
./install --all
```

After that, add to `~/.vimrc` the following lines:

```
...
Plugin 'Valloric/YouCompleteMe'
...

" ycm config
let g:ycm_global_ycm_extra_conf = '~/.vim/bundle/YouCompleteMe/third_party/ycmd/.ycm_extra_conf.py'
set completeopt=longest,menu
autocmd InsertLeave * if pumvisible() == 0|pclose|endif
inoremap <expr> <CR>       pumvisible() ? "\<C-y>" : "\<CR>"
inoremap <expr> <Down>     pumvisible() ? "\<C-n>" : "\<Down>"
inoremap <expr> <Up>       pumvisible() ? "\<C-p>" : "\<Up>"
let g:ycm_min_num_of_chars_for_completion=3
let g:ycm_seed_identifiers_with_syntax=1
"end

```

The documentation w.r.t. to the `.ycm_extra_conf.py` is shit and the location of it changes with diff versions. The location I found is diff from stackoverflow...
Right now this setting works for Cpp.
