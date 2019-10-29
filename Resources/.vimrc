" An example for a vimrc file.
"
" Maintainer:	Bram Moolenaar <Bram@vim.org>
" Last change:	2001 Jul 18
"
" To use it, copy it to
"     for Unix and OS/2:  ~/.vimrc
"	      for Amiga:  s:.vimrc
"  for MS-DOS and Win32:  $VIM\_vimrc
"	    for OpenVMS:  sys$login:.vimrc

" When started as "evim", evim.vim will already have done these settings.
if v:progname =~? "evim"
  finish
endif

" Use Vim settings, rather then Vi settings (much better!).
" This must be first, because it changes other options as a side effect.
set nocompatible

" allow backspacing over everything in insert mode
set backspace=indent,eol,start

set autoindent		" always set autoindenting on

" START ## MY SETUP - atle

" i want smartcase searching - i.e. if i write all letters with
" the same case, case doesnt matter, but it does if i do mix
" them
set ignorecase
set smartcase
" 4 SPACES INSTEAD OF TABS
set expandtab
set smarttab
set tabstop=4
set softtabstop=4
set shiftwidth=4

" SHOW LINE NUMBER
set number

" COMMAND SHELL IS CYGWIN-TCSH
"set shell=C:\\cygwin\\usr\\bin\\tcsh.exe DOES NOT WORK!!

" MAKE IT EASIER TO JUMP BETWEEN FILES (USING CTRL-J / CTRL-K)
map <C-J> <C-W>j<C-W>_
map <C-K> <C-W>k<C-W>_

" MAKES A MINIMAL 'WINDOW' REALLY SMALL...
set wmh=0

" JAVA WINHELP DOCUMENTATION FOLLOWS...
" RIGHT-CLICK
":amenu PopUp.JavaHelp :!start winhlp32 -k <cword> C:\Docs\Java\JavaWinhelp\JDK14.HLP <CR>
" SHORT CUT (NORMAL MODE)
"map <C-X> :!start winhlp32 -k <cword> C:\Docs\Java\JavaWinhelp\JDK14.HLP <CR>

" MAP FOR ctags. i use F3-F5 (cannot make the standard Ctrl-] work...
" F5 is the most useful!!
" use ctags -R . to create tags-file
se nocscopetag   " very important!! override what is in /etc/vimrc
" jump to other file. use Ctrl-t to get back
map <F3> :exec("tag ".expand("<cword>"))<CR>
" open in new tab
map <F4> :tab split<CR>:exec("tag ".expand("<cword>"))<CR>
" open in new vertical tab - more useful
map <F5> :vsp <CR>:exec("tag ".expand("<cword>"))<CR>
"map <C-]> :ta <cword> <CR> "!dir <cword> <CR> " DOES NOT WORK. WHY ??? (dir <cword> works...)
"nmap <c-r> :ta expand("<cword>")<CR>
"map <c-r> :ta expand("<cword>")<CR>
"map <F3> :ta <cword> <CR>
"map <F3> :ta expand("<cword>")<CR>
"map <C-]> <C-S>

"insert blank line is F1
:map <F1> o<Esc>

" show paranthesis mathes
se sm

" this is clipped from :help hex
"  Change that "*.bin" to whatever
"  comma-separated list of extension(s) you find yourself wanting to edit.
"  vim -b : edit binary using xxd-format!
augroup Binary
  au!
  au BufReadPre  *.fvp let &bin=1
  au BufReadPost *.fvp if &bin | %!xxd
  au BufReadPost *.fvp set ft=xxd | endif
  au BufWritePre *.fvp if &bin | %!xxd -r
  au BufWritePre *.fvp endif
  au BufWritePost *.fvp if &bin | %!xxd
  au BufWritePost *.fvp set nomod | endif
augroup END

" make sure pasting works fine in terminals (see :help paste)
" no, more useful to have indenting correct. -atle,
"set paste

" omnicompletion stuff
filetype plugin on
set ofu=syntaxcomplete#Complete
" maybe not so useful now that Ctrl-X Ctrl-O works...
let g:pydiction_location = '/private/agy/.vim/after/ftplugin/pydiction/complete-dict'

" for python debugging
:map <F2> opdb.set_trace()<Esc>



" END   ## MY SETUP -atle

if has("vms")
  set nobackup		" do not keep a backup file, use versions instead
else
  set backup		" keep a backup file
endif
set history=50		" keep 50 lines of command line history
set ruler		" show the cursor position all the time
set showcmd		" display incomplete commands
set incsearch		" do incremental searching

" For Win32 GUI: remove 't' flag from 'guioptions': no tearoff menu entries
" let &guioptions = substitute(&guioptions, "t", "", "g")

" Don't use Ex mode, use Q for formatting
map Q gq

" Make p in Visual mode replace the selected text with the "" register.
vnoremap p <Esc>:let current_reg = @"<CR>gvs<C-R>=current_reg<CR><Esc>

" This is an alternative that also works in block mode, but the deleted
" text is lost and it only works for putting the current register.
"vnoremap p "_dp

" Switch syntax highlighting on, when the terminal has colors
" Also switch on highlighting the last used search pattern.
if &t_Co > 2 || has("gui_running")
  syntax on
  set hlsearch
endif

" Only do this part when compiled with support for autocommands.
if has("autocmd")

  " Enable file type detection.
  " Use the default filetype settings, so that mail gets 'tw' set to 72,
  " 'cindent' is on in C files, etc.
  " Also load indent files, to automatically do language-dependent indenting.
  filetype plugin indent on

  " For all text files set 'textwidth' to 78 characters.
  " autocmd FileType text setlocal textwidth=78

  " When editing a file, always jump to the last known cursor position.
  " Don't do it when the position is invalid or when inside an event handler
  " (happens when dropping a file on gvim).
  autocmd BufReadPost *
    \ if line("'\"") > 0 && line("'\"") <= line("$") |
    \   exe "normal g`\"" |
    \ endif
  
  " make sure tabstop = 4 is respected for python
  autocmd FileType python setlocal shiftwidth=4 tabstop=4 softtabstop=4 expandtab

endif " has("autocmd")

"" auto-completion not so efficient...
"" want to have auto-completion of paranthesis'es etc.
" :inoremap ( ()<Esc>i
" :inoremap [ []<Esc>i
" :inoremap { {}<Esc>i
" :inoremap ' ''<Esc>i
" :inoremap " ""<Esc>i

" want to fold functions etc
set foldmethod=indent
set foldlevel=99
" Enable folding with the spacebar
"nnoremap <space> za

" aligns text (like comments after code) so they end at column 80
" <leader> means '\', so \+tab executes this
nnoremap <leader><tab> mc80A <esc>080lDgelD`cP
