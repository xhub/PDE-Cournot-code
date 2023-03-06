# SPDX-License-Identifier: MIT
using Logging, LoggingExtras, Crayons


function _my_fmt(level, _module, group, id, file, line)
  if level == Logging.Info
    return :blue, "", ""
  else
    return Logging.default_metafmt(level, _module, group, id, file, line)
  end
end

function get_str(cw::Crayons.CrayonWrapper)::String
  return join(get_str.(cw.v))
end

function get_str(s::String)::String
  return s
end

function get_logger(io::IO)

  simple_io_logger = 
  FormatLogger(io) do _io, args
    # remove ANSI escape characters
    # XXX: ugly thing
    println(_io, get_str(args.message))
#    println(_io, replace(args.message.v[1], r"[^\x20-\x7e]" => ""))
  end

  demux_logger = TeeLogger(
    simple_io_logger,
    ConsoleLogger(stdout; meta_formatter=_my_fmt)
  )

  return demux_logger
end
