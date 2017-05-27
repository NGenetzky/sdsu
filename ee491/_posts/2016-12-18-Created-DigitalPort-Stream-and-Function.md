---
title: Created DigitalPort, Stream and Function
layout: post
tags:
- gabriel
- particle
- matlab
source-id: 1OnW6E9injsD21ZcGL0Htflg3XweqlrTiVBDNszY1QXU
published: true
---
**Topic:** Create classes to encapsulate one or many digital pins.

**Start Time:** (12:00pm)	**Stop Time:** (04:30pm)	**Total Time Spent:** (4.5hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** I need to make combine Digital Input/Output classes.

* * *


**Class ****DigitalPin**

Combining DigitalInput & DigitalOutput

Good: Less classes. PinMode can be changed without destructing the class.

Bad: There is only a single vector (actually std::bitset) of Digital Pins. Therefore we are limited to get/set of 32 Digital Pins and the "set" will be useless on an Input.

**Class ****DigitalPort**

Store a vector of Digital Pins. Provide access to each of the get/set functions. Provides a easy constructor that will avoid all the "add( (DigitalPin) board_led)" sort of calls.

It is now much easier to create an app with certain Digital Pins exposed to the cloud api.

This is a 4 step process. I will show what I did to provide access to the digital pins on the Freenove smart car remote.

1. Create Pin and specify the direction (similar to pinMode or setting TRIS)

![image alt text](/sdsu/public/2MPVX1QHi24xa6XZM34BA_img_0.png)

![image alt text](/sdsu/public/2MPVX1QHi24xa6XZM34BA_img_1.png)

2. Create a Digital Pin.

![image alt text](/sdsu/public/2MPVX1QHi24xa6XZM34BA_img_2.png)

3. Create a Digital Port from the Digital Pins you want access to.

![image alt text](/sdsu/public/2MPVX1QHi24xa6XZM34BA_img_3.png)

4. Add the Digital Port to your App.

![image alt text](/sdsu/public/2MPVX1QHi24xa6XZM34BA_img_4.png)

**Topic:** Use Particle Device from Matlab

**Start Time:** (04:30pm)	**Stop Time:** (06:30pm)	**Total Time Spent:** (2.0hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** Need to use REST API from Matlab

* * *


[Earlier](https://docs.google.com/document/d/1MfsHqb2PtRfqkllonqriMbOFVn3qAxwobkmKheZnfAQ/edit?usp=sharing) we showed how to call these functions using the [Particle Command Line Interface](https://github.com/spark/particle-cli).

Now lets use it from **Matlab**!

Earlier I wrote a Matlab function called Particle. Given an ACCESS_TOKEN it can find all the devices associated and expose the functions and variables it has. It uses the REST API provided by particle.

Ideally it should have been a class rather than a function. Instead I use anonymous functions to make the object returned from the function behave as though it was a class. This is because I am very weak at Matlab Object Oriented Programming.

First call the Particle function to get anonymous functions to use "get" and “set”. The ACCESS_TOKEN for truncated for security reasons.

![image alt text](/sdsu/public/2MPVX1QHi24xa6XZM34BA_img_5.png)

Using set we can turn on all the LEDs. This call will have no effect on the Digital Pins configured as inputs.

	![image alt text](/sdsu/public/2MPVX1QHi24xa6XZM34BA_img_6.png)

Using get we can read the state of the Digital Pins. This will read the both pins configured as Inputs and Outputs. The switches are active low. The Digital Port configuration is shown on [1].

All the LEDs are OFF and all switches OFF:

![image alt text](/sdsu/public/2MPVX1QHi24xa6XZM34BA_img_7.png)

**Topic:** Use Particle functions string of characters.

**Start Time:** (06:00pm)	**Stop Time:** (11:30pm)	**Total Time Spent:** (5.5hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** Want to be able to parse function from stream of characters.

* * *


Next I want to provide this same functionality to a device that is connect to the microcontroller via Serial. I want to make it very generic where the commands come from but will assume that all commands will come as a stream of characters. All responses should should also be sent to a stream of characters.

**Class ****Stream**

Wraps a std::vector of characters.

Good: By using a vector I can easily obtain a char* that represents the data (elements of a vector are stored contiguously).

Bad: I have to regularly clear the buffer. The char* returned by data() will move if the container resizes.

Should be easy to use with the special "serialEvent()" function provided by arduino.

I will create two streams for my application. They are named std_in and std_out; these are named after stdin and stdout from C++ have similar purposes. The following are exposed to the cloud:

	Variables: stdin, stdout				Functions: cin

Now the task of interpreting command and data from the stdin, calling the appropriate functions with arguments and providing the response on stdout.

**Class ****Function**

Handle parsing Function from strings or vectors of characters.

I will choose to use bash's syntax for calling a function. That is: "$(function_name arguments)".

Constraints: Function name can't contain a space or ')’. Arguments can’t contain ‘)’.

I was able to successfully:

parse a function from stdin

Look up the function from a map

Execute the function

Place the result on stdout

	

	For instance, I read the switch values (with all switches off). 112 = 01110000

<table>
  <tr>
    <td>└─> particle call Logan cin "\$(get)"
6
└─> particle get Logan stdin
$(get)
└─> particle get Logan stdout
112</td>
  </tr>
</table>


Unfortunately the code for this is currently very ugly. It needs to be refactored and classes need to be created to make it easier.

