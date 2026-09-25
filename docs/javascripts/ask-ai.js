(() => {
  const root = document.getElementById("ask-ai");
  if (!root) return;

  const status = document.getElementById("ask-ai-status");
  const messages = document.getElementById("ask-ai-messages");
  const form = document.getElementById("ask-ai-form");
  const input = document.getElementById("ask-ai-input");
  const send = document.getElementById("ask-ai-send");
  const local = ["localhost", "127.0.0.1"].includes(location.hostname);
  const endpoint = "http://127.0.0.1:8765";
  const history = [];
  const languageAliases = {
    "c++": "cpp", cython: "python", f90: "fortran", f95: "fortran",
    fortran90: "fortran", py: "python", pyx: "python", sh: "bash",
    shell: "bash", yml: "yaml",
  };

  if (window.hljs) {
    hljs.configure({ languages: ["python", "fortran", "bash", "cpp", "c", "yaml", "javascript", "json", "xml"] });
  }

  function scrollToLatest() {
    messages.scrollTop = messages.scrollHeight;
  }

  function protectMath(markdown) {
    const expressions = [];
    const prefix = `ASKAIMATH${Math.random().toString(36).slice(2).toUpperCase()}TOKEN`;
    let source = "";
    let fence = null;
    let i = 0;
    while (i < markdown.length) {
      if (i === 0 || markdown[i - 1] === "\n") {
        const end = markdown.indexOf("\n", i);
        const line = markdown.slice(i, end === -1 ? undefined : end);
        const marker = /^ {0,3}(`{3,}|~{3,})/.exec(line);
        if (marker) {
          const character = marker[1][0];
          if (!fence) {
            fence = { character, length: marker[1].length };
          } else if (character === fence.character && marker[1].length >= fence.length &&
                     /^\s*$/.test(line.slice(marker[0].length))) {
            fence = null;
          }
          source += line;
          i += line.length;
          continue;
        }
      }
      if (fence) {
        source += markdown[i++];
        continue;
      }
      if (markdown[i] === "`") {
        const run = /^`+/.exec(markdown.slice(i))[0];
        const close = markdown.indexOf(run, i + run.length);
        if (close !== -1) {
          source += markdown.slice(i, close + run.length);
          i = close + run.length;
          continue;
        }
      }
      if (markdown[i] === "\\") {
        const close = markdown[i + 1] === "(" ? "\\)" : markdown[i + 1] === "[" ? "\\]" : null;
        if (close) {
          const end = markdown.indexOf(close, i + 2);
          if (end !== -1) {
            expressions.push(markdown.slice(i, end + 2));
            source += `${prefix}${expressions.length - 1}END`;
            i = end + 2;
            continue;
          }
        }
        source += markdown.slice(i, i + 2);
        i += 2;
        continue;
      }
      if (markdown[i] === "$") {
        const delimiter = markdown[i + 1] === "$" ? "$$" : "$";
        if (delimiter === "$$" || !/\s/.test(markdown[i + 1] || "")) {
          let end = i + delimiter.length;
          while ((end = markdown.indexOf(delimiter, end)) !== -1) {
            if (markdown[end - 1] !== "\\" &&
                (delimiter === "$$" || (markdown[end - 1] !== " " &&
                 !markdown.slice(i, end).includes("\n")))) break;
            end += delimiter.length;
          }
          if (end !== -1) {
            expressions.push(markdown.slice(i, end + delimiter.length));
            source += `${prefix}${expressions.length - 1}END`;
            i = end + delimiter.length;
            continue;
          }
        }
      }
      source += markdown[i++];
    }
    return { source, expressions, prefix };
  }

  function renderMarkdown(content, markdown) {
    if (!window.marked || !window.DOMPurify) {
      content.textContent = markdown;
      return;
    }
    const math = protectMath(markdown);
    const html = marked.parse(math.source, { gfm: true, breaks: true });
    content.innerHTML = DOMPurify.sanitize(html, {
      ALLOWED_TAGS: ["a", "blockquote", "br", "code", "del", "em", "h2", "h3", "h4", "hr",
        "li", "ol", "p", "pre", "strong", "table", "tbody", "td", "th", "thead", "tr", "ul"],
      ALLOWED_ATTR: ["href", "title", "class", "start", "align"],
    });
    if (math.expressions.length) {
      const pattern = new RegExp(`${math.prefix}(\\d+)END`, "g");
      const walker = document.createTreeWalker(content, NodeFilter.SHOW_TEXT);
      const nodes = [];
      while (walker.nextNode()) {
        if (walker.currentNode.textContent.includes(math.prefix)) nodes.push(walker.currentNode);
      }
      nodes.forEach((node) => {
        const fragment = document.createDocumentFragment();
        let from = 0;
        for (const match of node.textContent.matchAll(pattern)) {
          fragment.append(document.createTextNode(node.textContent.slice(from, match.index)));
          const expression = document.createElement("span");
          expression.className = "arithmatex";
          expression.textContent = math.expressions[Number(match[1])];
          fragment.append(expression);
          from = match.index + match[0].length;
        }
        fragment.append(document.createTextNode(node.textContent.slice(from)));
        node.replaceWith(fragment);
      });
    }
    content.querySelectorAll("a[href]").forEach((link) => {
      const url = new URL(link.href, location.href);
      if (!["http:", "https:"].includes(url.protocol)) {
        link.removeAttribute("href");
      } else if (url.origin !== location.origin) {
        link.target = "_blank";
        link.rel = "noopener noreferrer";
      }
    });
    content.querySelectorAll("pre > code").forEach((code) => {
      const language = [...code.classList].find((name) => name.startsWith("language-"))?.slice(9).toLowerCase() || "";
      const normalized = languageAliases[language] || language;
      if (window.hljs && hljs.getLanguage(normalized)) {
        code.className = `language-${normalized}`;
        hljs.highlightElement(code);
      }
      const block = document.createElement("div");
      block.className = "ask-ai__code-block";
      const label = document.createElement("div");
      label.className = "ask-ai__code-label";
      label.textContent = language || "code";
      block.append(label);
      code.parentElement.replaceWith(block);
      const pre = document.createElement("pre");
      pre.append(code);
      block.append(pre);
      const copy = document.createElement("button");
      copy.type = "button";
      copy.className = "ask-ai__copy ask-ai__copy-code";
      copy.textContent = "Copy code";
      copy.setAttribute("aria-label", "Copy code block");
      block.append(copy);
    });
  }

  function appendMessage(role, text = "") {
    const bubble = document.createElement("div");
    bubble.className = `ask-ai__message ask-ai__message--${role}`;
    if (role === "assistant") {
      const content = document.createElement("div");
      content.className = "ask-ai__message-content";
      bubble.append(content);
      renderMarkdown(content, text);
    } else {
      bubble.textContent = text;
    }
    messages.append(bubble);
    scrollToLatest();
    return bubble;
  }

  function setReply(bubble, text) {
    renderMarkdown(bubble.querySelector(".ask-ai__message-content"), text);
    scrollToLatest();
  }

  function setThinking(bubble, text) {
    let indicator = bubble.querySelector(".ask-ai__thinking");
    if (!text) {
      indicator?.remove();
      return;
    }
    if (!indicator) {
      indicator = document.createElement("div");
      indicator.className = "ask-ai__thinking";
      const circle = document.createElement("span");
      circle.className = "ask-ai__thinking-circle";
      circle.setAttribute("aria-hidden", "true");
      const label = document.createElement("span");
      indicator.append(circle, label);
      bubble.append(indicator);
    }
    indicator.lastElementChild.textContent = text;
    scrollToLatest();
  }

  async function typesetMath(bubble) {
    if (!window.MathJax) return;
    const content = bubble.querySelector(".ask-ai__message-content");
    try {
      if (MathJax.startup?.promise) await MathJax.startup.promise;
      if (!MathJax.typesetPromise) return;
      await MathJax.typesetPromise([content]);
      scrollToLatest();
    } catch (error) {
      console.warn("Could not render Ask AI math:", error);
    }
  }

  async function readStream(response, bubble) {
    if (!response.body) throw new Error("This browser does not support streamed replies.");
    const reader = response.body.getReader();
    const decoder = new TextDecoder();
    let buffer = "";
    let answer = "";
    let done = false;
    let frame = 0;
    const render = () => {
      frame = 0;
      setReply(bubble, answer);
    };
    const handleLine = (line) => {
      if (!line.trim()) return;
      const event = JSON.parse(line);
      if (event.type === "delta") {
        setThinking(bubble, "");
        answer += event.text;
        if (!frame) frame = requestAnimationFrame(render);
      } else if (event.type === "status") {
        setThinking(bubble, event.text);
      } else if (event.type === "done") {
        answer = event.answer;
        done = true;
        setThinking(bubble, "");
      } else if (event.type === "error") {
        throw new Error(event.message);
      }
    };
    try {
      while (true) {
        const { value, done: ended } = await reader.read();
        if (ended) break;
        buffer += decoder.decode(value, { stream: true });
        let newline;
        while ((newline = buffer.indexOf("\n")) !== -1) {
          handleLine(buffer.slice(0, newline));
          buffer = buffer.slice(newline + 1);
        }
      }
      buffer += decoder.decode();
      if (buffer.trim()) handleLine(buffer);
      if (!done) throw new Error("The reply stopped before it finished.");
      return answer;
    } finally {
      if (frame) cancelAnimationFrame(frame);
      setReply(bubble, answer);
      reader.releaseLock();
    }
  }

  async function checkHealth() {
    if (!local) {
      status.textContent = "Ask AI is available in the local documentation preview only.";
      input.disabled = true;
      send.disabled = true;
      return;
    }
    try {
      const response = await fetch(`${endpoint}/health`);
      if (!response.ok) throw new Error("Unavailable");
      const data = await response.json();
      const ready = data.ready && data.protocol === "ndjson-v1";
      if (data.protocol !== "ndjson-v1") {
        status.textContent = "Restart the Ask AI server to enable streaming replies.";
      } else if (!data.ready) {
        status.textContent = "Local server is running. Set OPENAI_API_KEY to enable chat.";
      } else if (ready) {
        const reasoning = data.reasoning_effort ? ` · ${data.reasoning_effort} reasoning` : "";
        status.textContent = `Local assistant ready · ${data.model}${reasoning} · ${data.commit}`;
      }
      status.dataset.ready = String(ready);
      input.disabled = !ready;
      send.disabled = !ready;
    } catch {
      status.textContent = "Start the local Ask AI server, then reload this page.";
      input.disabled = true;
      send.disabled = true;
    }
  }

  messages.addEventListener("click", async (event) => {
    const button = event.target.closest(".ask-ai__copy-code");
    if (!button) return;
    const text = button.closest(".ask-ai__code-block").querySelector("code").textContent;
    try {
      await navigator.clipboard.writeText(text);
      const previous = button.textContent;
      button.textContent = "Copied";
      setTimeout(() => { if (button.isConnected) button.textContent = previous; }, 1600);
    } catch {
      button.textContent = "Copy failed";
    }
  });

  form.addEventListener("submit", async (event) => {
    event.preventDefault();
    const question = input.value.trim();
    if (!question || send.disabled) return;
    appendMessage("user", question);
    input.value = "";
    send.disabled = true;
    input.disabled = true;
    const bubble = appendMessage("assistant");
    setThinking(bubble, "Thinking…");
    try {
      const response = await fetch(`${endpoint}/chat`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ message: question, history: history.slice(-8) }),
      });
      if (!response.ok) {
        const data = await response.json();
        throw new Error(data.detail || `Request failed (${response.status})`);
      }
      const answer = await readStream(response, bubble);
      history.push({ role: "user", content: question }, { role: "assistant", content: answer });
      await typesetMath(bubble);
    } catch (error) {
      setThinking(bubble, "");
      setReply(bubble, `Could not answer: ${error.message}`);
    } finally {
      send.disabled = false;
      input.disabled = false;
      input.focus();
    }
  });

  input.addEventListener("keydown", (event) => {
    if (event.key === "Enter" && !event.shiftKey && !event.isComposing) {
      event.preventDefault();
      form.requestSubmit();
    }
  });
  checkHealth();
})();
