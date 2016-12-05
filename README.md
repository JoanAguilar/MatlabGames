Two Simple Matlab Games
=======================

Overview
--------
Included are two simple videogames that are coded and run using Matlab. I coded these games  as a fun side-project with the intention to use them as teaching tool for an undergrad Matlab programming class. Unfortunately, I was not able to finish them on time. Both games use a limited number of Matlab functionalities and are less than 1000 lines long when combined. A second motivation was to test Matlab capabilities to run a videogame that included real time graphics.


Bang Bang!
----------
*Bang Bang!* is a 2-player aim game that puts each player in control of a cannon. The game progresses in turns and the goal is to hit the other player's cannon with a projectile before they hit you. The game includes realistic projectile physics (including wind).

![](https://github.com/JoanAguilar/Matlab-games/blob/master/images/Bangbang.png "Bang Bang! gameplay")

To play the game, run the file *Bangbang.m* using Matlab or Octave. You will first be asked to choose the terrain and wind conditions. From this point, players will take turns to shoot a projectile. The players get to select the initial angle and velocity of the projectile.

*Bang Bang!* was inspired by the 1990 Windows game with the same name. You can play the original [here](http://playdosgamesonline.com/bang-bang.html).


Lander
------
*Lander* is a 2D Moon landing simulation single player game. The goal consists of landing the Moon lander safely using a limited amount of fuel. In order to land safely, the vertical and horizontal velocities, the terrain angle, and the lander angle, all need to be less than a certain threshold. The user has control of the main thruster and the RCS thrusters (to control orientation). The game includes realistic physics with the lander modeled after the Apollo missions Moon lander.

![](https://github.com/JoanAguilar/Matlab-games/blob/master/images/Lander.png "Lander gameplay")

To play *Lander* run the file *Lander.m* using Matlab (recommended) or Octave. You will first be asked to select which controls you want to use and the level of difficulty. Use the controls you selected to land the Moon lander safely on the Moon.

*Lander* was inspired by the 1990 Windows game with the same name. You can find an Android version of the original [here](https://play.google.com/store/apps/details?id=com.pilot51.lander).


Disclaimer
----------
Feel free to use any of the code here for your own purposes. However, I will appreciate if you reference this repository when you do so :smiley:. If you want, you can also let me know what you are using it for, I will be glad to know, but this is optional.
