---
title: Communicate with photon via WiFi through cloud
layout: post
tags:
- matlab
- gabriel
- particle
- tinker
- thingspeak
source-id: 1h21xta-Wy08qsY5a6NRLQuCtsizAjbgwOSXdVOTogBk
published: true
---
**Topic:** Demo previous work

**Start Time:** (10:30pm)	**Stop Time:** (01:30pm)	**Total Time Spent:** (1.5hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** Need to demonstrate that I fulfilled half of the goal for the 2nd half of EE491.

* * *


# Communicate with Photon via WiFi through cloud

Tasks

1. Read a digital value

2. Write a digital value

3. Read an analog value

I will demonstrate this functionality in multiple ways.

1. Bash Script using Particle CLI.

2. Matlab using REST API

3. Tinker Android App

4. Graphs on Webinterface using Thingspeak

## Bash Script using Particle CLI

First I will first see what functions and variables are provided by my device.

<table>
  <tr>
    <td>+ particle list Logan
Logan [39001c001847353236343033] (Photon) is online
  Variables:
    help (string)
    stdin (string)
    stdout (string)
    data (string)
  Functions:
    int cin(String args) 
    int get(String args) 
    int set(String args) 
    int digitalread(String args) 
    int digitalwrite(String args) 
    int analogread(String args) 
    int analogwrite(String args) 
    int reg(String args) </td>
  </tr>
</table>


We will do the three tasks in order.

1. Read a digital value

    1. Get entire DigitalPort

<table>
  <tr>
    <td>+ particle call Logan get ' '
120</td>
  </tr>
</table>


    2. Read a DigitalPin from DigitalPort. 

        1. Note that D4 does not indicate D4 on Microcontroller. Instead each pin from (D0 to C8 correspond to a DigitalPin we added to the DigitalPort.

<table>
  <tr>
    <td>+ particle call Logan digitalread D4
1</td>
  </tr>
</table>


2. Write a digital value

    3. Write to the entire DigitalPort

<table>
  <tr>
    <td>+ particle call Logan set 97
0</td>
  </tr>
</table>


    4. Write to a single DigitalPin of the DigitalPort

<table>
  <tr>
    <td>+ particle call Logan digitalwrite D0=LOW
1</td>
  </tr>
</table>


3. Read an analog value

    5. Read using RegisterBank interface

<table>
  <tr>
    <td>+ particle call Logan reg 2
4095</td>
  </tr>
</table>


    6. Read using AnalogRead function

<table>
  <tr>
    <td>+ particle call Logan analogread A0
4095</td>
  </tr>
</table>


## Matlab using REST API

First I obtain a struct using Particle.m:

<table>
  <tr>
    <td>>> P = Particle(ACCESS_TOKEN)
P = 
    Comrad: [1x1 struct]
       Eli: [1x1 struct]
    Parker: [1x1 struct]
      Matt: [1x1 struct]
     Logan: [1x1 struct]</td>
  </tr>
</table>


Each of these structs contain anonymous functions that will interact with the REST API.

name.call_fx(args) 	-> P.call(name,fx,args)	 -> particle call name fx args

name.get_var() 	-> P.get(name,var)		 -> particle get name var

We will be using Logan, which is connected to a Freenove board (4 A Inputs, 3 D Out, 3 D Inputs)

<table>
  <tr>
    <td>>> P.Logan
ans = 
                 info: [1x1 struct]
                 name: 'Logan'
                   id: '39001c001847353236343033'
            connected: 1
                  url: 'https://api.particle.io/v1/devices/39001c001847353236343033/'
               has_fx: @(fx)any(strcmp(Y.(name).info.functions,fx))
              has_var: @(var)any(strcmp(fieldnames(Y.(name).info.variables),var))
             call_cin: @(args)P.call(Y.(name),fx,args)
             call_get: @(args)P.call(Y.(name),fx,args)
             call_set: @(args)P.call(Y.(name),fx,args)
     call_digitalread: @(args)P.call(Y.(name),fx,args)
    call_digitalwrite: @(args)P.call(Y.(name),fx,args)
      call_analogread: @(args)P.call(Y.(name),fx,args)
     call_analogwrite: @(args)P.call(Y.(name),fx,args)
             call_reg: @(args)P.call(Y.(name),fx,args)
             get_help: @(args)P.get(Y.(name),var)
            get_stdin: @(args)P.get(Y.(name),var)
           get_stdout: @(args)P.get(Y.(name),var)
             get_data: @(args)P.get(Y.(name),var)</td>
  </tr>
</table>


To understand how to use this device we will get the "help" variable.

<table>
  <tr>
    <td>>> help = P.Logan.get_help(); disp(help.result)
EE491 Particle Microcontroller
Variables:
help // How to use Application.
stdin // Characters in input buffer.
stdout // Characters in output buffer.
Data // String representing RegisterBank
Functions:
cin('string') -> stdin.append(string)
get('i') -> DigitalPort::get()
set('v') -> DigitalPort::set(v)
digitalread('LN') -> DigitalPort::set(N+L#)
digitalwrite('LN=D') -> DigitalPort::set(N+L#,D)
analogread('LN') -> Register::get(N+L#)
analogwrite('LN=A') -> Register::set(N+L#,A)
reg('i') -> Register::get(i);
reg('i=v') -> Register::set(i,v);
Notes:
// N=[0,8], B={HIGH,LOW}, A=[0,4095]
// L={A|B|C|D}, L#={D=0,A=8,B=16,C=24};
d0=DigitalPort( {board_led, led1, led2, led3, sw1, sw2, sw3} );
regs=RegisterBank( {t, d0, a0, a1, a2, a3} );</td>
  </tr>
</table>


We will do the three tasks in order.

4. Read a digital value

    7. Get entire DigitalPort

<table>
  <tr>
    <td>>> P.Logan.call_get(''); disp(ans.return_value);
    97
>> disp(de2bi(ans.return_value)); % [sw3, sw2, sw1, led3, led2, led1, board_led]
     1     0     0     0     0     1     1
</td>
  </tr>
</table>


    8. Read a DigitalPin from DigitalPort. 

        2. Note that D4 does not indicate D4 on Microcontroller. Instead each pin from (D0 to C8 correspond to a DigitalPin we added to the DigitalPort.

<table>
  <tr>
    <td>>> P.Logan.call_digitalread('D4'); disp(ans.return_value); % sw1
     0</td>
  </tr>
</table>


5. Write a digital value

    9. Write to the entire DigitalPort

<table>
  <tr>
    <td>>> P.Logan.call_set(num2str(bi2de([ 1 1 1 1 0 1 1 ])));; % [sw3, sw2, sw1, led3, led2, led1, board_led]</td>
  </tr>
</table>


    10. Write to a single DigitalPin of the DigitalPort

<table>
  <tr>
    <td>>> P.Logan.call_digitalwrite('D0=LOW');</td>
  </tr>
</table>


6. Read an analog value

    11. Read using RegisterBank interface

<table>
  <tr>
    <td>>> P.Logan.call_reg('2'); disp(ans.return_value); % POT1-A0
        4095</td>
  </tr>
</table>


    12. Read using AnalogRead function

<table>
  <tr>
    <td>>> P.Logan.call_analogread('A0'); disp(ans.return_value); % POT1-A0
        4095</td>
  </tr>
</table>


I wrote a Matlab function (ParticleDevice) that is designed to read all the variables, registers, and digital pins.

<table>
  <tr>
    <td>>> logan = ParticleDevice(P.Logan)
logan = 
    register: [8161716 110 4095 1346 941 2875]
       dport: [0 1 1 1 0 1 1]
      stdout: ''
       stdin: ''
        help: 'EE491 Particle Microcontroller
Variables:
help // How to use Application.
stdin // Characters in input buffer.
stdout // Charact...'
        data: '{"1":"0008164361","2":"0000000110","3":"4095","4":"1359","5":"1223","6":"2879",}'</td>
  </tr>
</table>


## Tinker Android App

The screenshot below shows the pins available on a Photon Microcontroller. The Particle team developed the app to communicate with firmware they designed called "Tinker". The firmware is designed to received commands in a certain format from the android app.

I developed Tinker.h, which allows me to define my own behavior for these commands. Each pin can have up to four commands: analogRead, analogWrite, digitalRead, digitalWrite.

Intuitively I connected the digital commands each with a DigitalPin. The analog commands work well with Register. A0-A3 are simply connected to Registers related to analog pins. A4 shows time since boot, and A5 and DAC are connected with the DigitalPort. 

<table>
  <tr>
    <td>



Nothing
DigitalPort
DigitalPort
millis()
JOYX-A3
JOYY-A2
POT2-A1
POT1-A0
</td>
    <td></td>
    <td>



Nothing
sw3
sw2
sw1
led3
led2
led1
board_led
</td>
  </tr>
</table>


## Graphs on Web Interface using Thingspeak

(Screenshots are from previous entry on 2017-12-24)

![image alt text]({{ site.url }}/public/KWVsE55nJKMAfqtQW6YgQ_img_0.png)

![image alt text]({{ site.url }}/public/KWVsE55nJKMAfqtQW6YgQ_img_1.png)

![image alt text]({{ site.url }}/public/KWVsE55nJKMAfqtQW6YgQ_img_2.png)

