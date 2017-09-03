# 极简主义
In virtualbox, first enable efi mode...    
use BIOS in vbox:   

```
parted /dev/sda
    mklabel msdos
    mkpart primary ext4 1M 100M
    set 1 boot on
    mkpart primary ext4 100M -1
    p
    q
```


```
ping -c 3 github.com
timedatectl set-ntp true
ls /sys/firmware/efi/efivars

fdisk -l
parted /dev/sdx
enter parted:
    mklabel gpt
    mkpart primary 1 512M (/boot)
    mkpart primary 512M 25G (/)
    mkpart primary 25G -1 (reserve)
    set 1 boot on
    p 
    q
leave parted.
mkfs.ext4 /dev/sda2
mkfs.vfat /dev/sda1

mount /dev/sda2 /mnt
mkdir -p /mnt/boot/efi
mount /dev/sda1 /mnt/boot/efi

nano /etc/pacman.d/mirrorlist
pacman -Syy

pacstrap -i /mnt base base-devel
genfstab -U -L /mnt >> /mnt/etc/fstab

arch-chroot /mnt /bin/bash
pacman -Syu zsh vim
zsh
ln -s /usr/share/zoneinfo/Asia/Shanghai /etc/localtime
hwclock --systohc --utc
alias ls='ls --color'
nano /etc/locale.gen
locale-gen
echo "LANG=en_US.UTF-8" > /etc/locale.conf
echo 'steven-mba' > /etc/hostname
nano /etc/hosts (replace myhostname)

mkinitcpio -p linux
passwd
useradd -m -g users -G wheel -s /bin/bash steven 
passwd steven
visudo

pacman -S grub-bios efibootmgr
grub-install --target=x86_64-efi --efi-directory=/boot/efi --bootloader-id=arch_grub --recheck --debug
grub-mkconfig -o /boot/grub/grub.cfg

pacman -Syu dialog wpa_supplicant
exit
umount -R /mnt
reboot
```

图形界面

```
lspci | grep -e VGA
xorg-server xorg-xrdb dbus xf86-video-vesa virtualbox-guest-utils

```

