---
title: Training materials
tags: TraintheTrainers, Talk
description: PARIS, ELIXIR-TtT 2020.
---

## Training materials: sharing and making re-use possible

<!-- Put the link to this slide here so people can follow -->
slide: https://hackmd.io/@Vic2020/rki40XkXL#/

---

**Learning outcome:** 
Be able to identify training materials that exist already, and develop a routine of sharing training materials with best practices

---

### When you prepare the Training material

#### You have two possibilities
- Starting from scratch
- Re-use and re-adapt existing Training material

---

### 70% of trainers in life sciences :heart: TeSS.

In order to prepare good material
Inspire you from other people material
- Check the licenses or ask permission
- Cite your sources
- Acknowledge (if the cases that you don't have references)

---

![][ELIXIR TeSS](https://tess.elixir-europe.org/)

---

### Usage FAIR Principles

---
![](https://i.imgur.com/VRLvkjU.png)

---

### Architecture of extension

---

![](https://i.imgur.com/ij69tPh.png)

---

## Content script

- Bind with each page
- Manipulate DOM
- Add event listeners
- Isolated JavaScript environment
  - It doesn't break things

---

# :fork_and_knife: 

---

<style>
code.blue {
  color: #337AB7 !important;
}
code.orange {
  color: #F7A004 !important;
}
</style>

- <code class="orange">onMessage('event')</code>: Register event listener
- <code class="blue">sendMessage('event')</code>: Trigger event

---

# :bulb: 

---

- Dead simple API
- Only cares about application logic

---

```typescript
import * as Channeru from 'channeru'

// setup channel in different page environment, once
const channel = Channeru.create()
```

---

```typescript
// in background script
const fakeLogin = async () => true

channel.answer('isLogin', async () => {
  return await fakeLogin()
})
```

<br>

```typescript
// in inject script
const isLogin = await channel.callBackground('isLogin')
console.log(isLogin) //-> true
```

---

# :100: :muscle: :tada:

---

### Wrap up

- Cross envornment commnication
- A small library to solve messaging pain
- TypeScript Rocks :tada: 

---

### Thank you! :sheep: 

You can find me on

- GitHub
- Twitter
- or email me
