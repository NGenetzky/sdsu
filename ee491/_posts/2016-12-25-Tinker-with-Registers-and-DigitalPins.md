---
title: Tinker with Registers and DigitalPins
layout: post
tags:
- gabriel
- particle
- tinker
source-id: 1_xL1B_F6Z7VRJOKgXOHjGL8w7VZDxXLMbA7NvdtJ1dk
published: true
---
**Topic:** Tinker with Registers and DigitalPins

**Start Time:** (12:00pm)	**Stop Time:** (01:30pm)	**Total Time Spent:** (1.5hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** 

* * *


On Christmas Eve I was in a rush to get to dinner. I shutdown my workspace before pushing my commits to github. I lost a few hours of work from yesterday. 

First I rewrote Register and RegisterBank from memory. I then split the "thingspeak" portion of RegisterBank into a class called FixedFields.

**Class ****FixedFields**

When constructing the object you must specify the number of characters for each field.

Can set values by calling a "set" function with an index of the field.

Has setup function to create json object that Particle's webhook understands.

Can get a string representation of it, or publish the string.

This could also be useful for parsing objects from Strings.

**Class ****Tinker**

	Designed to make it very easy to respond to commands from the [Tinker Android App](https://play.google.com/store/apps/details?id=io.particle.android.app&hl=en).

	I first generalized the 4 functions so that they can use the same function prototype.

	The class holds a vector of anonymous functions called "handlers".

When a command is sent from the app, each function calls the generalized function prototype. Each of the handlers will then be called with the values until request is handled.

WKP is not connected, DAC and A5 write to DigitalPort. A4 is ms since boot. A0-A3 read analog ports (two pots and joystick) . D0 is board LED, D1-D3 are LEDs. D4-D7 read the switches. 

Screenshot of Tinker App using the Tinker class

![image alt text](/sdsu/public/yveu3DGt0Pe859TwHqZclQ_img_0.png)

Freenove Smart Car Remote

![image alt text](/sdsu/public/yveu3DGt0Pe859TwHqZclQ_img_1.png)

