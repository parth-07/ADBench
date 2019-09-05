using System;
using Xunit;

namespace DotnetRunnerTests
{
    public class RunnerTests
    {
        [Fact]
        public void LibraryLoadTest()
        {
            var modulePath = "./MockTest.dll";
            var moduleLoader = new DotnetRunner.ModuleLoader(modulePath);
            Assert.True(moduleLoader.GetGMMTest() != null);
            Assert.True(moduleLoader.GetBATest() != null);
            Assert.True(moduleLoader.GetHandTest() != null);
            Assert.True(moduleLoader.GetLSTMTest() != null);

        }
    }
}
