===========================================
 Installation and running IPython notebook
===========================================

E-Cell4 docker for Ubuntu
===============================

Please install docker with http://docs.docker.io/en/latest/installation/ubuntulinux/

::

   sudo docker pull kozo2/ecell4-docker
   sudo docker run -p 8123:8888 -t kozo2/ecell4-docker
   # please open 0.0.0.0:8123 with your browser on Ubuntu host
   

E-Cell4 docker for Mac
===============================

Please install boot2docker and docker with http://docs.docker.io/en/latest/installation/mac/ .

And create and run following shell script.

::

   # vm must be powered off
   for i in {8000..8900}; do
    VBoxManage modifyvm "boot2docker-vm" --natpf1 "tcp-port$i,tcp,,$i,,$i";
    VBoxManage modifyvm "boot2docker-vm" --natpf1 "udp-port$i,udp,,$i,,$i";
   done

Next, please run following commands

::

   boot2docker up
   boot2docker ssh
   # now you are in boot2docker
   docker pull kozo2/ecell4-docker
   docker run -p 8123:8888 -t kozo2/ecell4-docker
   # please open 0.0.0.0:8123 with your browser on Mac
