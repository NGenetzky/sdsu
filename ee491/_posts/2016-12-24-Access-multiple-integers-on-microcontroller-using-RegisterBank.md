---
title: Access multiple integers on microcontroller using RegisterBank
layout: post
tags:
- particle
- thingspeak
- gabriel
source-id: 1Q_SC3GUFbKhhee7alBq42jYUvod6EE4ROQnGZFAnCVY
published: true
---
**Topic:** Access multiple integers on microcontroller using RegisterBank

**Start Time:** (11:00pm)	**Stop Time:** (04:00pm)	**Total Time Spent:** (5.0hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** Would like to provide easy access to variables on the device.

* * *


I wanted to be able to read the analog values from the micro controller. Instead of programming for this exact case I just said I wanted to read an integer from the microcontroller and that the value should be updated when I get it.

**Class ****Register**

Stores an int.  Provides get/set methods to get the cached value.

Stores an anonymous function for b

By default get/set will call the anonymous functions to update the cached value.	

**Class ****RegisterBank**

	Holds a vector of Registers. Provide get/set that accept index of register to act on.

	Hold vector of characters to hold string representation and provides Particle Variable.

	Stores index where each register should be printed.

	update() will call get new value for the register and update the chars for that register.

I used anonymous functions to allows the user to define how to "get" or “set” the integer value.

I defined these for current time (in ms), the DigitalPort mentioned earlier, and 4 analog values.

![image alt text]({{ site.url }}/public/vPOeHmqamX0AxJGBQ3l2eg_img_0.png)

Character representation of RegisterBank can be easily configured

	We would like to use it with a particle webhook and therefore we need a json object.

Current it is a json object with register index as key and register value as value.

Create a Thingspeak Channel

![image alt text]({{ site.url }}/public/vPOeHmqamX0AxJGBQ3l2eg_img_1.png)

Define a webhook that will send data to Thingspeak

![image alt text]({{ site.url }}/public/vPOeHmqamX0AxJGBQ3l2eg_img_2.png)

Publish an event to update the values on Thingspeak

name="thingspeak”, data=string representation of RegisterBank:

<table>
  <tr>
    <td>{"data":"{\"1\":\"0000577304\",\"2\":\"0112\",\"3\":\"2796\",\"4\":\"1135\",\"5\":\"1135\",\"6\":\"3054\"}","ttl":"60","published_at":"2016-12-24T21:41:49.295Z","coreid":"39001c001847353236343033","name":"thingspeak"}</td>
  </tr>
</table>


Visualize the values in Charts:

![image alt text]({{ site.url }}/public/vPOeHmqamX0AxJGBQ3l2eg_img_3.png)

![image alt text]({{ site.url }}/public/vPOeHmqamX0AxJGBQ3l2eg_img_4.png)

![image alt text]({{ site.url }}/public/vPOeHmqamX0AxJGBQ3l2eg_img_5.png)

