# Ask AI

Ask questions about Photochem's documentation, source code, and companion `photochem_clima_data` repository. The assistant can search tracked text files and inspect HDF5 dataset structure, but cannot run Photochem or execute commands.

<div class="ask-ai" id="ask-ai" aria-label="Photochem AI assistant">
  <div class="ask-ai__status" id="ask-ai-status" role="status">Checking local assistant…</div>
  <div class="ask-ai__messages" id="ask-ai-messages" role="log" aria-live="polite">
    <div class="ask-ai__message ask-ai__message--assistant">What would you like to know about Photochem?</div>
  </div>
  <form class="ask-ai__form" id="ask-ai-form">
    <label for="ask-ai-input">Your question</label>
    <textarea id="ask-ai-input" rows="3" maxlength="4000" placeholder="Ask about the code, API, or tutorials…" required></textarea>
    <div class="ask-ai__actions">
      <span>Enter to send · Shift+Enter for a new line</span>
      <button class="md-button md-button--primary" type="submit" id="ask-ai-send">Send</button>
    </div>
  </form>
</div>
