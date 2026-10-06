// lib/util.nf

// Nextflow passes `--flag false` on the command line as the String "false",
// which is truthy in Groovy. Every boolean param must be read through this.
def asBool(v) {
    if (v instanceof Boolean) return v
    if (v == null) return false
    def s = v.toString().trim().toLowerCase()
    if (s in ['true',  't', 'yes', 'y', '1']) return true
    if (s in ['false', 'f', 'no',  'n', '0', '']) return false
    throw new IllegalArgumentException("Cannot interpret '${v}' as a boolean (use true or false)")
}
