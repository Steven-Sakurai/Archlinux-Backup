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

# load root partition
mount /dev/sda2 /mnt
# create bootloader
mkdir -p /mnt/boot/efi
mount /dev/sda1 /mnt/boot/efi

# update mirrorlist
nano /etc/pacman.d/mirrorlist
pacman -Syy

# installation
pacstrap -i /mnt base base-devel
genfstab -U -L /mnt >> /mnt/etc/fstab

# enter the root partition
arch-chroot /mnt /bin/bash
# install the basics stuff
pacman -Syu zsh vim dialog wpa_supplicant
zsh

# set locale and time
ln -s /usr/share/zoneinfo/Asia/Shanghai /etc/localtime
hwclock --systohc --utc
nano /etc/locale.gen
locale-gen
echo "LANG=en_US.UTF-8" > /etc/locale.conf
echo 'steven-mba' > /etc/hostname
nano /etc/hosts (replace myhostname)

mkinitcpio -p linux

# set password and accounts
passwd
useradd -m -g users -G wheel -s /bin/bash steven 
passwd steven
visudo

# grub2
pacman -S grub-bios efibootmgr
grub-install --target=x86_64-efi --efi-directory=/boot/efi --bootloader-id=arch_grub --recheck --debug
grub-mkconfig -o /boot/grub/grub.cfg

# finished leave the partition
exit
umount -R /mnt
reboot
```

# graphics

```
lspci | grep -e VGA
xorg-server xorg-xrdb dbus xf86-video-vesa virtualbox-guest-utils

```

