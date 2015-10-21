syntax on
colorscheme evening
filetype plugin indent on

behave mswin
vnoremap <Tab> >
vnoremap <S-Tab> <
so $VIMRUNTIME/mswin.vim
nmap <A-LeftMouse> ms<LeftMouse><C-v>`so
imap <A-LeftMouse> <Esc><C-v>`^ms<Esc>gi<LeftMouse><C-o><C-v>`so
vmap <A-LeftDrag> <LeftDrag>
vmap <A-LeftMouse> <C-v><LeftMouse>msgv`s
vmap <S-LeftMouse> v<LeftMouse>msgv`s
set mouse=ra

set showmatch
set smarttab
set smartindent
set ai
set wrap
set tabstop=4
set shiftwidth=4
"set foldenable
"set foldmethod=syntax
autocmd BufWinLeave *.* mkview
autocmd BufWinEnter *.* silent loadview 


map! <F5> <esc>:wa<cr>:make<cr>
map! <F6> <esc>:wa<cr>:make run<cr>
map! <F7> <esc>:wa<cr>:make valgrind<cr>
map! <F8> <esc>:wa<cr>:make cgdb<cr>
map! <F2> <esc>:set invnumber<cr>
map! <F3> <esc>:wa<cr>:prev<cr>
map! <F4> <esc>:wa<cr>:next<cr>

map <F5> :wa<cr>:make<cr>
map <F6> :wa<cr>:make run<cr>
map <F7> :wa<cr>:make valgrind<cr>
map <F8> :wa<cr>:make cgdb<cr>
map <F2> :set invnumber<CR>
map <F3> :wa<cr>:prev<cr>
map <F4> :wa<cr>:next<cr>

nnoremap <Space> za

set numberwidth=1
highlight LineNr term=bold cterm=NONE ctermfg=DarkGrey ctermbg=NONE gui=NONE guifg=DarkGrey guibg=NONE
