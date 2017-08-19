set nocompatible              
set rtp+=~/.vim/bundle/Vundle.vim
call vundle#begin()

Plugin 'VundleVim/Vundle.vim'
Plugin 'scrooloose/nerdtree'
Plugin 'ervandew/supertab'
Plugin 'easymotion/vim-easymotion'
Plugin 'Shougo/neocomplete.vim'
Plugin 'Shougo/vimshell.vim'
Plugin 'Raimondi/delimitMate'
Plugin 'vim-airline/vim-airline'
Plugin 'vim-airline/vim-airline-themes'

call vundle#end()            
filetype plugin indent on    

set bg=light
set backspace=2
set number
set autoread
set wildmenu
set ruler
set smartcase
set hlsearch
set magic
set incsearch
set lazyredraw
set encoding=utf8
set expandtab
set smarttab
set shiftwidth=4
set tabstop=4
set mouse=a
set clipboard=unnamed
set whichwrap+=<,>,h,l,[,]
syntax enable
colo molokai 


"Manage nerdtree
"show nerdtree with C-n
map <C-n> :NERDTreeToggle<CR>
"quit nerdtree automatically
autocmd bufenter * if (winnr("$") == 1 && exists("b:NERDTree") && b:NERDTree.isTabTree()) | q | endif
"show hidden files
let NERDTreeShowHidden=1
"fi

"Manage neocomplete
" Disable AutoComplPop.
"let g:acp_enableAtStartup = 0
" Use neocomplete.
let g:neocomplete#enable_at_startup = 1
" Use smartcase.
let g:neocomplete#enable_smart_case = 1
" Set minimum syntax keyword length.
let g:neocomplete#sources#syntax#min_keyword_length = 3
"fi

