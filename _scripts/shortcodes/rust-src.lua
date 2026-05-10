-- _scripts/shortcodes/rust-src.lua
-- Usage in .qmd:
--   {{< rust-src masbayes bayesa.rs >}}
--   {{< rust-src masreml gblup/mod.rs#L42 >}}
-- Renders an inline link to the file in the package's GitHub repo.

local REPOS = {
  masbayes = "https://github.com/bowo1698/masbayes/blob/main/src/rust/src/",
  masreml  = "https://github.com/bowo1698/masreml/blob/main/src/rust/src/",
}

return {
  ["rust-src"] = function(args, kwargs, meta)
    if #args < 2 then
      return pandoc.RawInline("html",
        '<span style="color:red">[rust-src: needs &lt;package&gt; &lt;path&gt;]</span>')
    end
    local pkg  = pandoc.utils.stringify(args[1])
    local path = pandoc.utils.stringify(args[2])
    local base = REPOS[pkg]
    if not base then
      return pandoc.RawInline("html",
        '<span style="color:red">[rust-src: unknown package "' .. pkg .. '"]</span>')
    end
    local url   = base .. path
    local label = "🦀 " .. pkg .. "/" .. path
    return pandoc.RawInline("html",
      '<a class="rust-source-link" href="' .. url ..
      '" target="_blank" rel="noopener">' .. label .. '</a>')
  end,
}
