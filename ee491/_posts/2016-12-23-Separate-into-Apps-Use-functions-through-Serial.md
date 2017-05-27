---
title: Separate into Apps. Use functions through Serial.
layout: post
tags:
- gabriel
- particle
source-id: 1B2gTJZRzd_kH33cNRD4VxqU5aBDGEqdjVeUZIHC0vg0
published: true
---
**Topic:** Separate into Apps. Use functions through Serial.

**Start Time:** (10:30am)	**Stop Time:** (10:30pm)	**Total Time Spent:** (12hr)

**Author: **Nathan Genetzky		[Nathan@Genetzky.us](mailto:Nathan@Genetzky.us)		[github.com/NGenetzky](https://github.com/NGenetzky)

**Primary Issue:** App is going to quickly grow and make it harder to start new project.

* * *


**Class ****App**

The App class is currently tightly coupled with the different classes it uses.

I want to avoid using inheritance because it degrades performance and it would lead to multiple inheritance which is frowned upon.

I will utilize preprocessor directives (#define, #if) to conditionally add apps to App.

What constraints should be imposed on these Apps?

	They must be represented by a single member variable.

	The type and name of this variable should be accessible from macros.

	Should define SETUP macro, which should be called within App::setup()

	Should define APP macro, which should be called within the body of App.

AppDigitalPort.h

	Adds a DigitalPort member named dport to the App.

	Provides add(DigitalPort) and add(DigitalPin) functions to App.

	Provides Particle.Functions: get, set, digitalread, digitalwrite

This makes it very easy to add these Apps or remove them.

	Set APPNAME_EN = 0 to make the APP and SETUP be do nothing (be empty).

**Class ****Function**

	Previously this class merely parsed a name and args from a vector.

	Now the Function is callable without arguments and returns int.

Calling Function::bind(std::function<int< String>>) will result in this new function being called with the args provided by the member variable in Function.

How do you define which functions the microcontroller can perform?

	Create a map with the second member of the pair being a Particle.Function: 

	![image alt text](/sdsu/public/gePsyU1BctB8A3GJDUjPg_img_0.png)

The std::bind(...) allows class methods to be called by providing "this" as first argument.

Now I want to call these functions through a serial interface.

I can use RealTerm to send and monitor serial communication via the virtual COM port

	![image alt text](/sdsu/public/gePsyU1BctB8A3GJDUjPg_img_1.png)

I can use the SerialHandle class I wrote in Matlab to send commands:

	![image alt text](/sdsu/public/gePsyU1BctB8A3GJDUjPg_img_2.png)

Renamed Stream to File

	Reason: The primary benefit is the contiguous memory. I feel like File better describes this.

The Stream capability is really provided by "read_cursor" which keeps track of the position within the file. I would like to refactor this out of the File class.

My next goal:

I would like to be able to read the analog pins at faster than a frequency of once per second. I would also like to provide the digital port value. 

Would like to publish to thingspeak servers in json format.

Would like to provide this as a Particle.variable.

