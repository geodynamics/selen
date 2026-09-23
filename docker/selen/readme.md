This container hosts a built version of selen.

docker run -it --rm -v $HOME/selen:/home/selen_user/work geodynamics/selen

This command will start the selen docker image and give you terminal access. Any changes made in the /home/selen_user/work directory will be reflected on the host machine at home/selen.
