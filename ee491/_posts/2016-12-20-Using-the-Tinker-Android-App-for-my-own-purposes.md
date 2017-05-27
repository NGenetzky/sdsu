---
title: Using the Tinker Android App for my own purposes
layout: post
tags:
- tinker
- particle
- gabriel
source-id: 1FXeGcEZWNs8HTtK1urkot22k3xkJqtGJE9vZIqmOW_4
published: true
---
**Topic:** Using the Tinker Android App for my own purposes

**Start Time:** (12:00pm)	**Stop Time:** (05:00pm)	**Total Time Spent:** (5.0hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** Would like to control my Particle device from the android app.

* * *


**Class ****Identifier**

	Has a const static map of words that it understands. (DICTIONARY)

	Provides fast conversion from Identifier (or unsigned) to string (char *).

	Should be used as the key in maps of objects.

**Overloading tinker functions**

Earlier it was mentioned that Particle has a "tinker" application.

	If this [tinker firmware](https://github.com/spark/firmware/blob/latest/user/applications/tinker/application.cpp) is programed on the Particle device.

	Then you can use the [Particle android app](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&cad=rja&uact=8&ved=0ahUKEwj8lOaM8IrRAhWD7iYKHau-A6oQFggoMAI&url=https%3A%2F%2Fplay.google.com%2Fstore%2Fapps%2Fdetails%3Fid%3Dio.particle.android.app%26hl%3Den&usg=AFQjCNH5QtkyHDRXxKAjDBDGbeMRYdbDTA&sig2=766kTB7SlEcJtOzTGrWg3w) to control the pins via the following functions.

![image alt text](/sdsu/public/KWFmLYL8pJ3qX3n27yUNw_img_0.png)

What does the android app expect from the firmware

	Syntax for args in general is "PIN=VALUE" for write or “PIN” for read.

	For example writing to D7 from android app sends "D7=HIGH"

	D can be replaced by A, B, C for Analog port, Election port B and C respectively.

I provided a different implementation of the DigitalWrite and DigitalRead functions. 

	I designed it to accept the same command that the previous functions accepted.

Instead of reading or writing to the digital pins it will send a "get"(for digitalread) or “set (for digitalwrite) to the DigitalPort class defined [earlier](https://docs.google.com/document/d/1OnW6E9injsD21ZcGL0Htflg3XweqlrTiVBDNszY1QXU/edit).

I can now use the android app to read or write anything that is defined as DigitalPin.

