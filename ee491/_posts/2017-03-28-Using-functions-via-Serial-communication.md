---
title: Using functions via Serial communication
layout: post
tags:
- realterm
- particle
- gabriel
source-id: 13pAgph1S-gc1kM0-VsKtEvEKJByLkDyTGUL-VTns25E
published: true
---
**Topic:** Using functions via serial communication

**Start Time:** (12:00pm)	**Stop Time:** (10:30pm)	**Total Time Spent:** (10.5hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** Need to demonstrate that I fulfilled half of the goal for the 2nd half of EE491.

* * *


# Communicate with Photon via Serial

The goal is to do the following tasks

1. Read a digital value

2. Write a digital value

3. Read an analog value

First we must setup the PC for serial communication

Find the COM port for your device

![image alt text](/sdsu/public/IRVmiSEQM2z9Wb1BKTN1Wg_img_0.png)

Setting up RealTerm 

Display

![image alt text](/sdsu/public/IRVmiSEQM2z9Wb1BKTN1Wg_img_1.png)

Port:

![image alt text](/sdsu/public/IRVmiSEQM2z9Wb1BKTN1Wg_img_2.png)

Send:

![image alt text](/sdsu/public/IRVmiSEQM2z9Wb1BKTN1Wg_img_3.png)

We will do the three tasks in order. 

I will show the asci text sent over serial and then the next like is the response. Every message is terminated by a line feed (LF).

1. Read a digital value

    1. Get entire DigitalPort

<table>
  <tr>
    <td>$(get)
112</td>
  </tr>
</table>


    2. Read a DigitalPin from DigitalPort. 

        1. Note that D4 does not indicate D4 on Microcontroller. Instead each pin from (D0 to C8 correspond to a DigitalPin we added to the DigitalPort.

<table>
  <tr>
    <td>$(DR D4)
1</td>
  </tr>
</table>


2. Write a digital value

    3. Write to the entire DigitalPort

<table>
  <tr>
    <td>$(set 97)
0</td>
  </tr>
</table>


    4. Write to a single DigitalPin of the DigitalPort

<table>
  <tr>
    <td>$(DR D0=LOW)
1</td>
  </tr>
</table>


3. Read an analog value

    5. Read using RegisterBank interface

<table>
  <tr>
    <td>$(reg 2)
2801</td>
  </tr>
</table>


    6. Read using AnalogRead function

<table>
  <tr>
    <td>$(AR A0)
2802</td>
  </tr>
</table>


A screenshot of RealTerm (green=out, yellow=in).

![image alt text](/sdsu/public/IRVmiSEQM2z9Wb1BKTN1Wg_img_4.png)

