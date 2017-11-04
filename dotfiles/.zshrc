export PATH=$HOME/bin:/usr/local/bin:$PATH

export ZSH=/home/steven/.oh-my-zsh
ZSH_THEME="lambda"
# Uncomment the following line to enable command auto-correction.
#ENABLE_CORRECTION="true"
plugins=(git)
source $ZSH/oh-my-zsh.sh

export LANG=en_US.UTF-8

# Preferred editor for local and remote sessions
if [[ -n $SSH_CONNECTION ]]; then
   export EDITOR='vim'
else
   export EDITOR='vim'
fi

# Compilation flags
export ARCHFLAGS="-arch x86_64"

# ssh
export SSH_KEY_PATH="~/.ssh/rsa_id"

# Set personal aliases, overriding those provided by oh-my-zsh libs,
alias c='clear'
alias d='cd ~/Desktop'
alias pacman='sudo pacman'
alias bright-adjust='echo 5 | sudo tee /sys/class/backlight/acpi_video0/brightness'
alias solarized-dark='cp ~/solarized-master/xfce4-terminal/dark/terminalrc ~/.config/xfce4/terminal/terminalrc' 
alias killsogou='killall fcitx && killall sogou-qimpanel && killall sogou-qimpanel-watchdog'
alias start-sogou='fcitx && sogou-qimpanel'
alias Git='git add -A && git commit -m "From Home, Arch" && git push orgin master'

