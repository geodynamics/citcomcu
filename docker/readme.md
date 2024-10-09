This container hosts a built version of citcomcu.

docker run -it --rm -v $HOME/citcomcu:/home/citcomcu_user/work geodynamics/citcomcu

This command will start the citcomcu docker image and give you terminal access. Any changes made in the /home/citcomcu_user/work directory will be reflected on the host machine at home/citcomcu.
