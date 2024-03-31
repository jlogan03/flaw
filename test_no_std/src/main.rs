//! Building this module successfully guarantees that the library is no-std compatible

#![no_std]
#![no_main]

use core::panic::PanicInfo;

use flaw;

#[panic_handler]
fn panic(_info: &PanicInfo) -> ! {
    // We can't print, so there's not much to do here
    loop {}
}

#[no_mangle]
pub fn _start() -> ! {
    let mut filter = flaw::butter1(0.01);
    let _ = filter.unwrap().update(1.0);

    loop {} // We don't actually run this, just compile it
}
