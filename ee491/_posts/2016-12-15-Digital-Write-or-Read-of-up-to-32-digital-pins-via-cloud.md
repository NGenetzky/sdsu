---
title: Digital Write or Read of up to 32 digital pins via cloud
layout: post
tags:
- particle
- tinker
- gabriel
source-id: 1MfsHqb2PtRfqkllonqriMbOFVn3qAxwobkmKheZnfAQ
published: true
---
**Topic:** Digital Write or Read of up to 32 digital pins via cloud 

**Start Time:** (06:00pm)	**Stop Time:** (12:00pm)	**Total Time Spent:** (6.0hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** I need to write a value to a digital pin from the cloud (or read from).

* * *


What does Particle provide to help with this?

[Tinker](https://github.com/spark/firmware/blob/60f703bd4165a292f89124a7e3755e4eb37b1a03/user/applications/tinker/application.cpp) provides 4 cloud functions firmware for reading and writing analog and digital pins.

<table>
  <tr>
    <td>Particle.function("digitalread", tinkerDigitalRead);
Particle.function("digitalwrite", tinkerDigitalWrite);
Particle.function("analogread", tinkerAnalogRead);
Particle.function("analogwrite", tinkerAnalogWrite);</td>
  </tr>
</table>


They have an [Android App](https://play.google.com/store/apps/details?id=io.particle.android.app&hl=en) that can utilize these function control the board.

I believe I could create functions with the same name that expects the same syntax. I would then be able to use their android app to do almost anything I want.

My code base had expanded dramatically. I am going to start from a blank slate for this new portion of the study. I intend to create very object oriented code.

Classes I created today:

Pin - Stores a pin number. Could be used to make sure only one instance of each pin exists. Also could provide functions that typically accept a pin number.

DigitalInput - Has a Pin. Sets pin to input on setup, and then provides functions to read.

DigitalOutput - Has a Pin. Sets pin to output on setup, and then provides functions to write.

App - Provides a global data structure. Can provide Particle Functions or Variables. Should standardizes what should be included in a firmware application.

I created a simple program that controls 4 LEDs and 3 Switches.

	I can read or write to this hardware from firmware.

	Or I can read or write to this hardware from the cloud.

This is an excerpt of "particle list" which lists all devices and the Functions/Variables it provides.

![image alt text]({{ site.url }}/public/vQMWDgU1kUsBw2O9FIQ_img_0.png)

I can write the value on up to 32 Digital Outputs via this command.

![image alt text]({{ site.url }}/public/vQMWDgU1kUsBw2O9FIQ_img_1.png)

I can read the value on up to 32 Digital Inputs via this command.

![image alt text]({{ site.url }}/public/vQMWDgU1kUsBw2O9FIQ_img_2.png)

