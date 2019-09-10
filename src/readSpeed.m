function [speed] = readSpeed(prhPath)

prh = load(prhPath);
speed = prh.speed.JJ;

end