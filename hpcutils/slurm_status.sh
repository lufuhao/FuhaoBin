#!/bin/bash

systemctl status slurmctld && echo
systemctl status slurmd.service && echo
sinfo && echo
